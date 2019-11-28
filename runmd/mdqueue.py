#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import os
import multiprocessing
import multiprocessing.queues


class MdQueue:

    def __init__(self, name, queue, gpus, function, resources=""):
        self.__name = name
        self.__queue = queue
        self.__gpus = gpus
        self.__function = function
        self.__procs = []
        self.__work_dir = os.path.abspath(os.curdir)

    def start(self):
        lock_filename = "lock." + self.__name
        if os.path.isfile(lock_filename):
            raise Exception("Could not capture lock '"+lock_filename+"'\n\nRemove lock via\n\n    rm "+lock_filename + "\n")
        else:
            lock_file = open(lock_filename, "w")
            lock_file.write("This lock prevents from running several copies simulteneuosly.")
            lock_file.close()

        q = multiprocessing.queues.SimpleQueue()

        for item in self.__queue:
            q.put(item)
        function = self.__function
        work_dir = self.__work_dir
        class P(multiprocessing.Process):
            def __init__(self, gpu_id):
                super(P, self).__init__()
                self.gpu = gpu_id

            def run(self):
                while not q.empty():
                    item = q.get()
                    if item is not None:
                        os.chdir(work_dir)
                        os.environ["CUDA_VISIBLE_DEVICES"] = str(self.gpu)
                        function(item)

        self.__procs = [P(i) for i in self.__gpus]

        for p in self.__procs:
            p.start()

    def join(self):
        for p in self.__procs:
            p.join()
        lock_filename = "lock." + self.__name
        if os.path.isfile(lock_filename):
            os.remove(lock_filename)


