#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import os
import marshal
import time
from subprocess import Popen, PIPE, STDOUT
import subprocess
import sys
import shutil
import re


class MdQueue:

    def __init__(self, name, tasks, function, resources="nodes=1:ppn=1", refresh_time=0.5, mail_when_done=False, mail_on_error=True):
        self.__name = name
        self.__tasks = tasks
        self.__function = function
        self.__resources = resources
        self.__work_dir = os.path.abspath(os.curdir)

        home_dir = os.path.expanduser("~")
        self.__pickles_dir = home_dir+"/.runmd_tmp/"+self.__name+"_pickles"
        self.__i = 0
        self.__pbs_ids = []
        self.__pickles_dir_created = False
        self.__sleep_time = refresh_time
        self.__mail_when_done = mail_when_done
        self.__mail_on_error = mail_on_error
        print refresh_time

    def join(self):
        dots = 0
        print self.__sleep_time
        while True:

            for pbs_id in self.__pbs_ids[:]:
                self.__check_pbs_id(pbs_id)

            if len(self.__pbs_ids) == 0:
                print "All tasks are done."
                break
            nrunning_now = len(os.listdir(self.__pickles_dir))
            CURSOR_UP_ONE = '\x1b[1A'
            ERASE_LINE = '\x1b[2K'
            print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)
            print "Running %d tasks. No tasks left in queue. Waiting running tasks for termination%s" % (nrunning_now, "."*dots)
            dots = (dots + 1) % 4
            time.sleep(self.__sleep_time)
        os.rmdir(self.__pickles_dir)
        self.__send_done_report()

    def __submit_task(self):
        if self.__i >= len(self.__tasks):
            return ""

        script_template = """#!/usr/bin/python
import os
import sys

work_dir = "%s"
filename="%s"

sys.path.append(work_dir)
#print os.path.abspath(os.curdir)

import marshal, types, traceback
try:
    os.chdir(work_dir)
    code, args = marshal.load(open(filename,"rb"))
#    print code, args
    func = types.FunctionType(code, globals(), "func_name")
    func(args)
except Exception as e:
    traceback.print_exc()
    sys.stderr.write(e.message)
    raise e;
finally:
    os.remove(filename)

        """
        f = open(self.__pickles_dir + "/" + str(self.__i), "wb")
        marshal.dump((self.__function.func_code, self.__tasks[self.__i]), f)
        f.close()

        script = script_template % (self.__work_dir, self.__pickles_dir + "/" + str(self.__i))
        p = Popen(['qsub',
                   '-N', self.__name +"_"+ str(self.__i) ,
                   '-l', self.__resources,
                   '-o', "stdoe/"+self.__name+"_"+str(self.__i)+".stdout",
                   '-e', "stdoe/"+self.__name+"_"+str(self.__i)+".stderr"
                   ],
                  stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        output = p.communicate(input=script)[0]
        p.wait()
        if p.returncode != 0:
            raise UnboundLocalError("Could not submit task, qsub exit code = %d\n%s"%(p.returncode,output));
        self.__i += 1
        self.__pbs_ids.append(output.strip())
        return output

    def __get_pbs_info(self, pbs_id):
        try:
            qstat = subprocess.check_output(["qstat", "-f", pbs_id])
            info = dict([map(lambda x: x.strip(),entry.split(" = ")) for entry in re.findall("    [\w]+ = [^\n]*", re.sub("\n\t", "", qstat))])
            return info
        except Exception:
            return {}

    def __check_pbs_id(self, pbs_id):
        try:
            info = self.__get_pbs_info(pbs_id)
            if info["job_state"] == "C":  # completed
                self.__pbs_ids.remove(pbs_id)
                exit_status = info["exit_status"]
                if exit_status != "0":
                    print >> sys.stderr, "Job `%s` exited with non-zero status %s.\n" % (pbs_id, exit_status)
                    sys.stderr.flush()
                    self.__send_error_report(pbs_id)
        except KeyError:
            pass

    def process_via_nproc(self, nproc):
        if os.path.isdir(self.__pickles_dir):
            print >> sys.stderr, ("Remove '"+self.__pickles_dir+"' to process. Abort.")
            return
        else:
            os.mkdir(self.__pickles_dir)
            self.__pickles_dir_created = True

        if not os.path.isdir("stdoe"):
            os.mkdir("stdoe")

        dots = 0
        print "Queueing started."
        while True:
            for pbs_id in self.__pbs_ids[:]:
                self.__check_pbs_id(pbs_id)

            if self.__i == len(self.__tasks):
                break

            nrunning_now = len(self.__pbs_ids)
            if nrunning_now < nproc:
                dots = 0
                self.__submit_task()
            else:
                CURSOR_UP_ONE = '\x1b[1A'
                ERASE_LINE = '\x1b[2K'
                print(CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE)
                print "Running %d tasks. Rest %d waits in queue%s" % (nrunning_now, len(self.__tasks)-self.__i, "."*dots)
                dots = (dots + 1) % 4
            time.sleep(self.__sleep_time)

        self.join()

    def interrupt(self):
        import subprocess
        try:
            for pbs_id in self.__pbs_ids:
                try:
                    dot = pbs_id.find(".")
                    if dot != -1:
                        idx = pbs_id[:dot]
                    else:
                        idx = pbs_id
                    subprocess.call(["qdel", idx])
                finally:
                    pass
        finally:
            if self.__pickles_dir_created:
                shutil.rmtree(self.__pickles_dir, ignore_errors=True)

    def __send_error_report(self, pbs_id):
        import postman

        if not self.__mail_on_error:
            return

        try:
            info=self.__get_pbs_info(pbs_id)
            error_path = info["Error_Path"][info["Error_Path"].find("/"):]
            output_path = info["Output_Path"][info["Error_Path"].find("/"):]

            cerr = self.__tail(error_path)
            cout = self.__tail(output_path)

            subject = "Fatal error in queue '%s' " % self.__name
            text = self.__disclaimer_noreply()
            text += "Job %s exited with non-zero status %s\n\n" % (pbs_id, info["exit_status"])
            text += "tail -n 20 stdout:\n" + 35*"-" + "\n" + cout+35*"-" + "\n\n"
            text += "tail -n 20 stderr:\n" + 35*"-" + "\n" + cerr+35*"-" + "\n"
            text += self.__disclaimer_noreply()
            # print "Subject: %s" % subject
            # print "Message: %s" % text
            # print
            postman.send_mail(subject, text, files=[error_path, output_path])
        except Exception:
            print >> sys.stderr, "Could not send mail for some reason. \n"
            sys.stderr.flush()

    def __send_done_report(self, ):
        import postman
        if not self.__mail_when_done:
            return
        subject = "Queue '%s' is done" % self.__name
        text = self.__disclaimer_noreply()
        text += "Queue '%s' is successfully finished. Well done!\n" % self.__name
        text += self.__disclaimer_noreply()
        postman.send_mail(subject, text)

    def __disclaimer_noreply(self):
        return "\n*** This is an automatically generated email, please do not reply ***\n\n"

    def __tail(self, filename, lines=20):
        try:
            f = open(filename, "r")
            total_lines_wanted = lines

            BLOCK_SIZE = 1024
            f.seek(0, 2)
            block_end_byte = f.tell()
            lines_to_go = total_lines_wanted
            block_number = -1
            blocks = [] # blocks of size BLOCK_SIZE, in reverse order starting
                        # from the end of the file
            while lines_to_go > 0 and block_end_byte > 0:
                if (block_end_byte - BLOCK_SIZE > 0):
                    # read the last block we haven't yet read
                    f.seek(block_number*BLOCK_SIZE, 2)
                    blocks.append(f.read(BLOCK_SIZE))
                else:
                    # file too small, start from begining
                    f.seek(0,0)
                    # only read what was not read
                    blocks.append(f.read(block_end_byte))
                lines_found = blocks[-1].count('\n')
                lines_to_go -= lines_found
                block_end_byte -= BLOCK_SIZE
                block_number -= 1
            all_read_text = ''.join(reversed(blocks))
            return '\n'.join(all_read_text.splitlines()[-total_lines_wanted:])
        except:
            return "Could not read `%s`" % filename