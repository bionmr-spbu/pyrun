#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import re
import sys
import subprocess
import time
import socket
import glob
import pickle

from pytrj import *
from pdblib.num import *

from runmd import trjtool
import shutil


class AmberTleapRc:
    def __init__(self):
        self.__shell_thickness = None
        self.__structure = None
        self.__output_dir = None
        self.__basename = None
        self.__water_model = "SPCBOX"
        self.__leaprc = None
        self.__ss_bonds = []
        self.__commands_before_solvatation = []
        self.__commands_after_solvatation = []
        self.__custom_solvatation_script = None

    def get_cmd(self):
        if not os.path.isdir(self.__output_dir) or\
           self.__structure is None or \
           self.__basename is None or \
           self.__shell_thickness is None:
            raise Exception("Could not create tleap command with paramters" +
                            ", ".join(map(str, [self.__shell_thickness,
                                                self.__structure,
                                                self.__output_dir,
                                                self.__basename])))
        else:
            self.__prepare_tleaprc(self.__output_dir)
            return ["tleap", "-s", "-f", self.__output_dir+"/tleap.rc"]

    def set_leaprc(self, leaprc_path):
        if not os.path.isfile(leaprc_path):
            raise Exception("No such file'"+leaprc_path+"'")
        else:
            self.__leaprc = os.path.abspath(leaprc_path)

    def _set_output_dir(self, output_dir):
        if not os.path.isdir(output_dir):
            raise Exception("No such directory '"+output_dir+"'")
        else:
            self.__output_dir = os.path.abspath(output_dir)

    def get_output_dir(self):
        return self.__output_dir

    def set_structure(self, mol):
        self.__structure = mol

    def set_basename(self, basename):
        self.__basename = str(basename)

    def set_shell_thickness(self, shell_thickness):
        self.__shell_thickness = float(shell_thickness)

    def set_water_model(self,water_model):
        self.__water_model = water_model

    def add_command(self, command):  # DEPRECIATED, use add_command_before_solvatation instead
        """
        :param command: tleap command
        :type command: str
        """
        self.add_command_before_solvatation(command)

    def add_command_before_solvatation(self, command):
        """
        :param command: tleap command
        :type command: str
        """
        self.__commands_before_solvatation.append(command)

    def add_command_after_solvatation(self, command):
        """
        :param command: tleap command
        :type command: str
        """
        self.__commands_after_solvatation.append(command)

    def get_water_model(self):
        return self.__water_model

    def set_custom_solvatation_script(self, script):
        """
        :type script: str
        """
        self.__custom_solvatation_script =  script

    def __prepare_tleaprc(self, work_dir):
        if not os.path.isdir(work_dir) or\
                self.__basename is None or \
                self.__shell_thickness is None:
            raise Exception("Could not create tleap command with paramters" +
                            ", ".join(map(str, [self.__shell_thickness,
                                                work_dir,
                                                self.__basename])))
        else:
            tleaprc = open(work_dir+"/tleap.rc", 'w')
            basename = self.__output_dir+"/"+self.__basename
            pdb_filename = basename+".auxiliary.pdb"
            self.__structure.write(pdb_filename)
            tleaprc.writelines(open(self.__leaprc).readlines())

            for command in self.__commands_before_solvatation:
                tleaprc.write(command + "\n")

            if self.__custom_solvatation_script is None:
                tleaprc.write("wbox = loadpdb "+pdb_filename+"\n")
                tleaprc.write("solvateoct wbox "+self.__water_model+" " + str(self.__shell_thickness)+"\n")
            else:
                tleaprc.write(self.__custom_solvatation_script)

            for command in self.__commands_after_solvatation:
                tleaprc.write(command + "\n")

            for bond in self.__ss_bonds:
                tleaprc.write("bond wbox."+str(bond[0])+".SG wbox."+str(bond[1])+".SG" + "\n")

            tleaprc.write("""
addions wbox Na+ 0
saveamberparm wbox """+basename+""".prmtop """+basename+""".inpcrd
savepdb wbox """+basename+""".pdb
quit
        """)
            tleaprc.close()

    def get_resulting_topology_filename(self):
        return self.__get_resulting_base_filename()+".prmtop"

    def get_resulting_pdb_filename(self):
        return self.__get_resulting_base_filename()+".pdb"

    def get_initial_structure(self):
        return self.__structure

    def get_resulting_inpcrd_filename(self):
        return self.__get_resulting_base_filename()+".inpcrd"

    def __get_resulting_base_filename(self):
        if not os.path.isdir(self.__output_dir) or\
           self.__basename is None:
            raise Exception("Could not evaluate base filename with paramters" +
                            ", ".join(map(str, [self.__output_dir,
                                                self.__basename])))
        else:
            return self.__output_dir+"/"+self.__basename

    def unset_ss_bond(self, residue_serial_1, residue_serial_2):
        pair = (min(residue_serial_1, residue_serial_2), max(residue_serial_1, residue_serial_2))
        if pair in self.__ss_bonds:
            self.__ss_bonds.remove(pair)

    def unset_ss_bonds(self):
        self.__ss_bonds = []

    def get_ss_bonds(self):
        """
        Returns residue numbers of bonded cyxteins
        :return: list of residue numbers
        :rtype: list of int
        """
        return copy(self.__ss_bonds)


    def set_ss_bond(self, residue_serial_1, residue_serial_2):
        """
        :type residue_serial_1: int
        :type residue_serial_2: int
        """
        pair = (min(residue_serial_1, residue_serial_2), max(residue_serial_1, residue_serial_2))
        self.__ss_bonds.append(pair)


class AmberTopologyModificator:

    def __init__(self):
        self.__original_topology = None
        self.__command = ""



    def _set_original_topology(self, filename):
        self.__original_topology = filename

    def get_type_of(self, mask):
        import StringIO
        output = StringIO.StringIO("")

        from ParmedTools.parmed_cmd import ParmedCmd
        from ParmedTools.parmlist import ParmList
        from ParmedTools.exceptions import (ParmError, InterpreterError)

        amber_prmtop = ParmList()
        amber_prmtop.add_parm(self.__original_topology+".mod")

        command = StringIO.StringIO("printLJTypes " + mask)

        parmed_commands = ParmedCmd(amber_prmtop, stdin=command, stdout=output)
        parmed_commands.use_rawinput = 0
        # parmed_commands.interpreter = opt.interpreter
        parmed_commands.prompt = ''
        # Loop through all of the commands
        try:
            parmed_commands.cmdloop()
        except InterpreterError, err:
            sys.exit('%s: %s' % (type(err).__name__, err))
        except ParmError:
            # This has already been caught and printed. If it was re-raised, then
            # that means we wanted to exit
            sys.exit(1)

        # DANGEROUS PARSING  ^_^
        lines = output.getvalue().split("\n")
        import re
        return int(re.split("\s+",lines[3].strip())[-1])



    def add_command(self, command):
        """
        :type command: str
        """
        self.__command += command +"\n"

    def change_angle(self,
                     atom_1_sel,
                     atom_2_sel,
                     atom_3_sel,
                     force_constant,
                     equilibrium_angle):
        self.__command += "setAngle "+" ".join(map(str,
                                                   (atom_1_sel,
                                                    atom_2_sel,
                                                    atom_3_sel,
                                                    force_constant,
                                                    equilibrium_angle)
                                                   )) + "\n"

    def change_bond(self,
                    atom_1_sel,
                    atom_2_sel,
                    force_constant,
                    equilibrium_distance):
        self.__command += "setBond "+" ".join(map(str,
                                                  (atom_1_sel,
                                                   atom_2_sel,
                                                   force_constant,
                                                   equilibrium_distance)
                                                  )) + "\n"

    def change_vdw(self,
                   atom_1_sel,
                   atom_2_sel,
                   radius,
                   well_depth):
        self.__command += "changeLJPair "+" ".join(map(str,
                                                       (atom_1_sel,
                                                        atom_2_sel,
                                                        radius,
                                                        well_depth)
                                                       )) + "\n"

    def delete_dihedral(self,
                        atom_1_sel,
                        atom_2_sel,
                        atom_3_sel,
                        atom_4_sel):
        self.__command += "deleteDihedral "+" ".join(map(str,
                                                         (atom_1_sel,
                                                          atom_2_sel,
                                                          atom_3_sel,
                                                          atom_4_sel)
                                                         )) + "\n"\


    def write_to(self, output_filename):
        if os.path.exists(output_filename):
            if os.path.samefile(output_filename, self.__original_topology):
                raise Exception("Attempt to overwrite "+output_filename+"'")

        from ParmedTools.parmed_cmd import ParmedCmd
        from ParmedTools.parmlist import ParmList
        from ParmedTools.exceptions import (ParmError, InterpreterError)

        amber_prmtop = ParmList()
        amber_prmtop.add_parm(self.__original_topology)

        command = self.__form_parmed_command(output_filename)

        parmed_commands = ParmedCmd(amber_prmtop, stdin=command, stdout=sys.stdout)
        parmed_commands.use_rawinput = 0
        # parmed_commands.interpreter = opt.interpreter
        parmed_commands.prompt = ''
        # Loop through all of the commands
        try:
            parmed_commands.cmdloop()
        except InterpreterError, err:
            sys.exit('%s: %s' % (type(err).__name__, err))
        except ParmError:
            # This has already been caught and printed. If it was re-raised, then
            # that means we wanted to exit
            sys.exit(1)

    def __form_parmed_command(self, output_filename):
        import StringIO
        command = self.__command
        command += "outparm " + " " + output_filename + "\n"
        command += "quit\n"
        return StringIO.StringIO(command)


class AmberInputFile:
    def __init__(self):
        self.__parameters = {
            "imin": None,
            "irest": None,
            "ntx": None,
            "ntb": None,
            "iwrap": None,
            "ntt": None,
            "tempi": None,
            "ntp": None,
            "pres0": None,
            "taup": None,
            "cut": None,
            "ntr": None,
            "ntc": None,
            "ntf": None,
            "nstlim": None,
            "dt": None,
            "ntpr": None,
            "ntwx": None,
            "ntwr": None,
            "ioutfm": None,
            "nmropt": None
        }
        self.__help = {}
        self.__noe_restraints_filename = None
        self._noe_restraints = {}
        self._noe_angle_restraints = {}
        self._noe_dihedral_restraints = {}
        self.__pinned_residues = []
        self._pinned_residues_with_filter = []

    def unset_noe_restraints(self):
        self._noe_restraints.clear()
        self.__update_nmropt()

    def unset_noe_restraint(self, atom_id1, atom_id2):
        atom_id1, atom_id2 = min(atom_id1, atom_id2), max(atom_id1, atom_id2)
        self._noe_restraints.pop((atom_id1, atom_id2), None)
        self.__update_nmropt()

    def unset_noe_angle_restraint(self, atom_id1, atom_id2, atom_id3):
        atom_id1, atom_id3 = min(atom_id1, atom_id3), max(atom_id1, atom_id3)
        self._noe_angle_restraints.pop((atom_id1,atom_id2, atom_id3), None)
        self.__update_nmropt()

    def unset_noe_angle_restraints(self):
        self._noe_angle_restraints.clear()
        self.__update_nmropt()

    def set_noe_restraint(self, atom_id1, atom_id2, r1, r2, r3, r4, k2, k3, comment=""):
        atom_id1, atom_id2 = min(atom_id1, atom_id2), max(atom_id1, atom_id2)
        self._noe_restraints[(atom_id1, atom_id2)] = (float(r1),
                                                      float(r2),
                                                      float(r3),
                                                      float(r4),
                                                      float(k2),
                                                      float(k3),
                                                      str(comment))
        self.__update_nmropt()

    def set_noe_angle_restraint(self, atom_id1, atom_id2, atom_id3, r1, r2, r3, r4, k2, k3, comment=""):
        atom_id1, atom_id3 = min(atom_id1, atom_id3), max(atom_id1, atom_id3)
        self._noe_angle_restraints[(atom_id1, atom_id2, atom_id3)] = (float(r1),
                                                      float(r2),
                                                      float(r3),
                                                      float(r4),
                                                      float(k2),
                                                      float(k3),
                                                      str(comment))
        self.__update_nmropt()

    def set_noe_dihedral_restraint(self, atom_id1, atom_id2, atom_id3, atom_id4, r1, r2, r3, r4, k2, k3, comment=""):
        atom_id1, atom_id3 = min(atom_id1, atom_id3), max(atom_id1, atom_id3)
        self._noe_dihedral_restraints[(atom_id1, atom_id2, atom_id3, atom_id4)] = (float(r1),
                                                      float(r2),
                                                      float(r3),
                                                      float(r4),
                                                      float(k2),
                                                      float(k3),
                                                      str(comment))
        self.__update_nmropt()

    def write(self, output_filename):
        output = open(output_filename, "w")
        self.__set_noe_restraints_filename(os.path.splitext(output_filename)[0]+".restraints")
        output.write(" &cntrl\n")
        for parameter in sorted(self.__parameters):
            value = self.__parameters[parameter]
            if not (value is None):
                output.write("    %-10s = %10s, " % (parameter, str(value)))
                if (parameter, value) in self.__help:
                    output.write("   ! " + self.__help[(parameter, value)])
                output.write("\n")
        output.write("/\n")

        if len(self._noe_restraints) > 0 or len(self._noe_angle_restraints) > 0 or ( hasattr(self, '_noe_dihedral_restraints') and  len(self._noe_dihedral_restraints)):
             output.write(" &wt type='END' / \n")

        for pin in self.__pinned_residues:
            output.write("Pinned residues:\n")
            output.write("%lf \nRES %d %d\n" % pin)
            output.write("END\n")

        if hasattr( self, '_pinned_residues_with_filter' ):
            for pin in self._pinned_residues_with_filter:
                force_constant, filters_strings, residues_ranges = pin;
                output.write("Pinned residues with filter:\n")
                output.write("%lf\n"%force_constant)
                output.write("FIND\n")
                for fltr in filters_strings:
                    output.write("%s\n"%fltr)
                output.write("SEARCH\n")
                for res_range in residues_ranges:
                    output.write("RES %d %d \n"%res_range)
            if len(self._pinned_residues_with_filter)>0:
                output.write("END\n")

        if len(self._noe_restraints) > 0 or len(self._noe_angle_restraints) > 0 or ( hasattr(self, '_noe_dihedral_restraints') and  len(self._noe_dihedral_restraints)):
            self.__write_noe_restraints()
            self.__write_noe_angle_restraints(append=True)
            if hasattr(self, '_noe_dihedral_restraints')  : self.__write_noe_dihedral_restraints(append=True)

            output.write("DISANG="+self.__noe_restraints_filename+"\n")

        output.write("END\n")
        output.close()

    def set(self, key, value, help_string=None):
        self.__parameters[key] = value
        if help_string is not None:
            self.__help[(key, value)] = str(help_string)
        return self

    def get(self, key):
        result = str(self.__parameters[key])
        try:
            result = int(result)
        except ValueError:
            result = float(result)
        finally:
            return result

    def set_pinned_residues(self, force_constant, first, last):
        self.__pinned_residues.append((float(force_constant), int(first), int(last)))

    def set_pinned_residues_with_filter(self, force_constant, filters_strings, list_of_res_ranges):
        """
        :type force_constant: float or int
        :type filters_strings: list of str
        :type list_of_res_ranges: list of tuple
        """
        if not hasattr( self, '_pinned_residues_with_filter' ):
            self._pinned_residues_with_filter = []
        self._pinned_residues_with_filter.append((float(force_constant), filters_strings, list_of_res_ranges))

    def __set_noe_restraints_filename(self, restraints_filename):
        self.__noe_restraints_filename = os.path.abspath(restraints_filename)

    def unset_pinned_residues(self):
        self.__pinned_residues = []

    def __write_noe_restraints(self,append=False):
        mode = "a" if append else "w"
        output = open(self.__noe_restraints_filename, mode)
        if output is None:
            raise Exception("Could not open file for writing '%s'" % self.__noe_restraints_filename)
        for pair in self._noe_restraints:
            noe = self._noe_restraints[pair]
            params = noe[-1:] + pair + noe[:-1]
            output.write("""
 &rst  ! %s
  ixpk= 0, nxpk= 0, iat=%d, %d, r1= %f, r2= %f, r3= %f, r4= %f,
      rk2=%f, rk3=%f, ir6=1, ialtd=0,
 &end
 """ % params)
        output.close()

    def __write_noe_angle_restraints(self, append=False):
        mode = "a" if append else "w"
        output = open(self.__noe_restraints_filename, mode)
        if output is None:
            raise Exception("Could not open file for writing '%s'" % self.__noe_restraints_filename)
        for pair in self._noe_angle_restraints:
            noe = self._noe_angle_restraints[pair]
            params = noe[-1:] + pair + noe[:-1]
            output.write("""
 &rst  ! %s
  ixpk= 0, nxpk= 0, iat=%d, %d, %d, r1= %f, r2= %f, r3= %f, r4= %f,
      rk2=%f, rk3=%f, ir6=1, ialtd=0,
 &end
 """ % params)
        output.close()

    def __write_noe_dihedral_restraints(self, append=False):
        mode = "a" if append else "w"
        output = open(self.__noe_restraints_filename, mode)
        if output is None:
            raise Exception("Could not open file for writing '%s'" % self.__noe_restraints_filename)
        for pair in self._noe_dihedral_restraints:
            noe = self._noe_dihedral_restraints[pair]
            params = noe[-1:] + pair + noe[:-1]
            output.write("""
 &rst  ! %s
  ixpk= 0, nxpk= 0, iat=%d, %d, %d, %d, r1= %f, r2= %f, r3= %f, r4= %f,
      rk2=%f, rk3=%f, ir6=1, ialtd=0,
 &end
 """ % params)
        output.close()

    def __update_nmropt(self):
        if len(self._noe_restraints) == 0 and len(self._noe_angle_restraints)==0 and  not (hasattr(self, '_noe_dihedral_restraints') and len(self._noe_dihedral_restraints)!=0):
            self.set("nmropt", 0, "No nmr-type analysis will be done")
        else:
            self.set("nmropt", 1, "NMR restraints and weight changes will be read")


class MD:

    FATAL = -1
    ERROR = FATAL+1
    WARNING = ERROR+1
    INFO = WARNING+1
    TRACE = INFO+1
    DEBUG = TRACE+1

    __str_level = {FATAL:   " FATAL ",
                   ERROR:   " ERROR ",
                   WARNING: "WARNING",
                   INFO:    " INFO  ",
                   TRACE:   " TRACE ",
                   DEBUG:   " DEBUG "}

    __build_dir_name = "1_build"
    __minimization_dir_name = "2_minimization"
    __heat_dir_name = "3_heat"
    __equilibration_dir_name = "4_equilibration"
    __run_dir_name = "5_run"
    __pattern = "%05d"
    __dump_filename = "auto_dump.pickle"

    @staticmethod
    def load_state(filename):
        """
        Returns MD stated, stored in pickle file
        :param filename: MD pickle
        :return: MD loaded from pickle
        :rtype: MD
        """
        f = open(filename, "rb")
        if f:
            result = pickle.load(f)
            result._set_hostname()
            return result
        else:
            raise OSError("Could not open '%s' for reading." % filename)

    def __init__(self, name="__no_name__"):
        self.is_verbose = False
        self.log_level = None
        self.keep_netcdf = False
        self.keep_restart = True

        self.__mol = Mol()
        self.__work_dir = None
        self.__log_filename = None
        self.__topology_filename = None

        self.__id = name
        self.__hostname = socket.gethostname()
        self.__current_step = 0
        self.__reference_structure = None

        self.tleaprc = AmberTleapRc()
        self.topology_mod = AmberTopologyModificator()

        self.min1_parameters = AmberInputFile()
        self.min2_parameters = AmberInputFile()
        self.heat_parameters = AmberInputFile()
        self.equil_parameters = AmberInputFile()
        self.run_parameters = AmberInputFile()

        self.__restart_filename = None
        self.__pmemd_executable = None


    def build(self):
        self.log("Building parameters...")
        if not os.path.samefile(os.curdir, self.__work_dir):
            self.log("Work dir is unset", MD.FATAL)

        if not os.path.isdir(MD.__build_dir_name):
            self.__make_dir(MD.__build_dir_name)

        self.log("Writing tleap parameters...", MD.DEBUG)

        self.tleaprc._set_output_dir(self.get_build_dir())
        tleap_cmd = self.tleaprc.get_cmd()
        self.__run_cmd(tleap_cmd)

        self.topology_mod._set_original_topology(self.tleaprc.get_resulting_topology_filename())
        self.__topology_filename = os.path.abspath(self.tleaprc.get_resulting_topology_filename()+".mod")
        self.topology_mod.write_to(self.__topology_filename)
        self.set_restart_file(self.tleaprc.get_resulting_inpcrd_filename())
        self.__mol = Mol(self.tleaprc.get_resulting_pdb_filename())
        self.log("Building done.")

    def minimize(self):
        self.log("Minimizing...")
        if not os.path.samefile(os.curdir, self.__work_dir):
            self.log("Work dir is unset", MD.FATAL)

        min_dir = self.get_minimization_dir()

        if not os.path.isdir(min_dir):
            self.__make_dir(min_dir)

        self.min1_parameters.write(min_dir+"/"+"min_1.in")
        self.min2_parameters.write(min_dir+"/"+"min_2.in")

        if self.__reference_structure is not None:
            ref_struct = self.__reference_structure
        else:
            ref_struct = self.__restart_filename

        min_1_cmd = self.__pmemd_executable + [
            "-O",
            "-i",   min_dir + "/" + "min_1" + ".in",
            "-o",   min_dir + "/" + "min_1" + ".out",
            "-p",   self.__topology_filename,
            "-c",   self.__restart_filename,
            "-ref", ref_struct,
            "-r",   min_dir + "/" + "min_1" + ".rst",
            "-inf", min_dir + "/" + "min_1" + ".mdinfo",
            "-l",   min_dir + "/" + "min_1" + ".log"
        ]

        min_2_cmd = self.__pmemd_executable + [
            "-O",
            "-i",   min_dir + "/" + "min_2" + ".in",
            "-o",   min_dir + "/" + "min_2" + ".out",
            "-p",   self.__topology_filename,
            "-ref", ref_struct,
            "-c",   min_dir + "/" + "min_1" + ".rst",
            "-r",   min_dir + "/" + "min_2" + ".rst",
            "-inf", min_dir + "/" + "min_2" + ".mdinfo",
            "-l",   min_dir + "/" + "min_2" + ".log"
        ]

        self.__run_cmd(min_1_cmd)
        self.__run_cmd(min_2_cmd)

        self.set_restart_file(min_dir + "/" + "min_2" + ".rst")

        self.log("Minimizing done.")

    def heat(self):
        self.log("Heating...")
        if not os.path.samefile(os.curdir, self.__work_dir):
            self.log("Work dir is unset", MD.FATAL)

        heat_dir = self.get_heat_dir()

        if not os.path.isdir(heat_dir):
            self.__make_dir(heat_dir)

        self.heat_parameters.write(heat_dir+"/"+"heat.in")

        if self.__reference_structure is not None:
            ref_struct = self.__reference_structure
        else:
            ref_struct = self.__restart_filename

        heat_cmd = self.__pmemd_executable + [
            "-O",
            "-i",   heat_dir + "/" + "heat" + ".in",
            "-o",   heat_dir + "/" + "heat" + ".out",
            "-p",   self.__topology_filename,
            "-c",   self.__restart_filename,
            "-ref", ref_struct,
            "-x",   heat_dir + "/" + "heat" + ".nc",
            "-r",   heat_dir + "/" + "heat" + ".rst",
            "-inf", heat_dir + "/" + "heat" + ".mdinfo",
            "-l",   heat_dir + "/" + "heat" + ".log"
        ]

        self.__run_cmd(heat_cmd)
        self.set_restart_file(heat_dir + "/" + "heat" + ".rst")
        self.log("Heating done.")

    def equilibrate(self):
        self.log("Equilibration...")
        if not os.path.samefile(os.curdir, self.__work_dir):
            self.log("Work dir is unset", MD.FATAL)

        equil_dir = self.get_equilibration_dir()

        if not os.path.isdir(equil_dir):
            self.__make_dir(equil_dir)

        self.equil_parameters.write(equil_dir+"/"+"equil.in")

        if self.__reference_structure is not None:
            ref_struct = self.__reference_structure
        else:
            ref_struct = self.__restart_filename

        heat_cmd = self.__pmemd_executable + [
            "-O",
            "-i",   equil_dir + "/" + "equil" + ".in",
            "-o",   equil_dir + "/" + "equil" + ".out",
            "-p",   self.__topology_filename,
            "-c",   self.__restart_filename,
            "-ref", ref_struct,
            "-x",   equil_dir + "/" + "equil" + ".nc",
            "-r",   equil_dir + "/" + "equil" + ".rst",
            "-inf", equil_dir + "/" + "equil" + ".mdinfo",
            "-l",   equil_dir + "/" + "equil" + ".log"
        ]

        self.__run_cmd(heat_cmd)
        self.set_restart_file(equil_dir + "/" + "equil" + ".rst")
        self.log("Equilibration done.")

    def do_md_step(self):
        performed_step = self.__current_step + 1
        self.log("Running step "+str(performed_step)+"...")
        if not os.path.samefile(os.curdir, self.__work_dir):
            self.log("Work dir is unset", MD.FATAL)

        run_dir = self.get_run_dir()

        if not os.path.isdir(run_dir):
            self.__make_dir(run_dir)


        pattern = MD.__pattern
        suffix = pattern % performed_step

        # self.run_parameters.__set_noe_restraints_filename(run_dir + "/" + "run" + suffix + ".restraints")
        self.run_parameters.write(run_dir + "/" + "run" + suffix + ".in")

        if self.__reference_structure is not None:
            ref_struct = self.__reference_structure
        else:
            ref_struct = self.__restart_filename

        if isinstance(self.__pmemd_executable,str):
            self.__pmemd_executable = [self.__pmemd_executable]
 
        run_cmd = self.__pmemd_executable + [
            "-O",
            "-i",   run_dir + "/" + "run" + suffix + ".in",
            "-o",   run_dir + "/" + "run" + suffix + ".out",
            "-p",   self.__topology_filename,
            "-c",   self.__restart_filename,
            "-ref", ref_struct,
            "-x",   run_dir + "/" + "run" + suffix + ".nc",
            "-r",   run_dir + "/" + "run" + suffix + ".rst",
            "-inf", run_dir + "/" + "run" + suffix + ".mdinfo",
            "-l",   run_dir + "/" + "run" + suffix + ".log"
        ]

        self.__run_cmd(run_cmd)
        self.__run_trjtool_on(run_dir + "/" + "run" + pattern + ".nc", performed_step)

        if not self.keep_netcdf:
            os.remove(run_dir + "/" + "run" + suffix + ".nc")

        if not self.keep_restart:
            os.remove(self.__restart_filename)

        self.set_restart_file(run_dir + "/" + "run" + suffix + ".rst")
        self.__current_step = performed_step


        self.log("Step "+str(performed_step)+" done.")

        self.save_state(filename=MD.__dump_filename)


    def set_restart_file(self, restart_filename):
        if not os.path.isfile(restart_filename):
            self.log("No such file '%s' " % restart_filename, MD.FATAL)
        else:
            self.__restart_filename = os.path.abspath(restart_filename)

    def set_step_as_restart_file(self, step_number=-1):
        """
        :param step_number:
        :type step_number: int
        :return:
        """
        if step_number < 0:
            step_number = self.get_current_step() + step_number + 1

        self.set_restart_file(self.get_run_dir() + "/" + "run" + MD.__pattern % step_number + ".rst")

    def set_pmemd(self, pmemd):
        if isinstance(pmemd,(str,unicode)):
            self.__pmemd_executable = [pmemd]
        else:
            self.__pmemd_executable = pmemd

    def get_build_dir(self):
        return os.path.abspath(self.__work_dir + "/" + MD.__build_dir_name)

    def get_structure(self):
        return self.__mol

    def get_minimization_dir(self):
        return os.path.abspath(self.__work_dir + "/" + MD.__minimization_dir_name)

    def get_heat_dir(self):
        return os.path.abspath(self.__work_dir + "/" + MD.__heat_dir_name)

    def get_equilibration_dir(self):
        return os.path.abspath(self.__work_dir + "/" + MD.__equilibration_dir_name)

    def get_run_dir(self):
        return os.path.abspath(self.__work_dir + "/" + MD.__run_dir_name)

    def set_reference_structure(self, reference_filename):
        self.__reference_structure = reference_filename

    def set_work_dir(self, work_dir_path):
        self.log("Changing work directory")
        if not os.path.isdir(work_dir_path):
            self.__make_dir(work_dir_path)
        self.__work_dir = os.path.abspath(work_dir_path)
        os.chdir(work_dir_path)
        self.log("Work directory changed to '"+work_dir_path+"'")

    def set_log_filename(self, log_filename):
        """ Sets log filename in WORKING DIRECTORY. Path is stripped. """
        self.log("Changing log filename to '"+log_filename+"'")
        log_path, log_basename = os.path.split(os.path.abspath(log_filename))
        if len(log_path) != 0:
            f = open(log_filename, 'a')
            f.close()
            self.log("Log filename should not contain '/' ", MD.WARNING)
            self.__log_filename = os.path.basename(log_filename)

    def set_current_step(self, value):
        self.__current_step = value

    def get_current_step(self):
        return self.__current_step


    def get_frame_fast(self, step, frame):
        if step < 0:
            step = self.get_current_step() + step+1
        if frame < 0:
            frame = (self.run_parameters.get("nstlim")+self.run_parameters.get("ntwx") -1)/self.run_parameters.get("ntwx") + frame
        suffix = self.__pattern % step
        dat_filename = self.get_run_dir() + "/" + "run" + suffix + ".dat"

        def read_frame(frame):
            finp = open(dat_filename, 'rb')
            nitem,ndim,dtype = unpack('III', finp.read(3*4))
            info=[]
            size=20
            offset1 = finp.tell() + nitem*(4+size)

            data = np.memmap(finp,
                                  dtype='float32',
                                  mode='r',
                                  offset=offset1+(frame*len(self.get_structure().getats())) *3*4,
                                  shape=(nitem,ndim))
            putmat(self.__mol, data, skipnull=True)
            finp.close()

        read_frame(frame)

        return self.__mol

    def get_frame(self, step, frame):
        if step < 0:
            step = self.get_current_step() + step+1
        if frame < 0:
            frame = (self.run_parameters.get("nstlim")+self.run_parameters.get("ntwx") -1)/self.run_parameters.get("ntwx") + frame
        suffix = self.__pattern % step
        dat_filename = self.get_run_dir() + "/" + "run" + suffix + ".dat"
        bf = Binfile(dat_filename)
        if len(bf.data) <= frame:
            self.log("Too big frame number (%d >= %d) in file '%s' abort." % (frame, len(bf.data), dat_filename), MD.FATAL)
        coor = array(bf.data[frame])
        if coor.shape[0] != len(self.__mol.getats()):
            self.log('the matrix does not match atoms %d != %d'%(coor.shape[0], len(self.__mol.getats())), MD.FATAL)
        putmat(self.__mol, coor, skipnull=True)
        return self.__mol

    def __run_trjtool_on(self, pattern, number):
        extract_dict = {
            "reference_structure": self.tleaprc.get_resulting_pdb_filename(),
            "netcdf_pattern": pattern,
            "netcdf_number": number,
            "frame_stride": "1",
            "frame_first_number": "",
            "frame_last_number": "",
            "segment_ids": "",
            "residue_ids": "1-%d" % len(self.__mol.getreses()),
            "residue_names": "",
            "atoms": "all",
            "output_filename": os.path.splitext(pattern % number)[0]+".dat",
            "is_ascii": "false",
        }

        trjtool_content = trjtool.get_trjtool_input(extract_dict)

        trjinp_filename = "trjtool" + str(os.getpid()) + ".inp"
        trjinp = open(trjinp_filename, "w")
        trjinp.write(trjtool_content)
        trjinp.close()

        trjtool_cmd = ["trjtool", trjinp_filename]
        self.__run_cmd(trjtool_cmd)

        os.remove(trjinp_filename)

    def log(self, message, level=INFO):
        """
        :param message: message
        :type message: str
        :param level: level of message
        """
        if level > MD.DEBUG:
            level = MD.DEBUG
        if level < MD.FATAL:
            level = MD.FATAL
        if level <= self.log_level:
            hostname = "[" + self.__hostname + "]"
            cuda_devices= "[" + os.getenv("CUDA_VISIBLE_DEVICES","?") + "]"
            md_id = "[" + self.__id + "]"
            date_time = time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime())
            level_str = "["+MD.__str_level[level]+"]"
            log_file = open(self.__work_dir+"/"+self.__log_filename, "a")
            log_string = hostname + cuda_devices + md_id + date_time + level_str + ": " + message+"\n"
            log_file.write(log_string)
            if self.is_verbose:
                sys.stdout.write(log_string)
                sys.stdout.flush()
            log_file.close()
        if level == MD.FATAL:
            self.__on_fatal()

    def __run_cmd(self, cmd, abort_on_failure=True):
        """ Returns tuple (exit_code, command_output) """
        self.log("Running " + cmd[0] + "...", MD.TRACE)
        self.log("with parameters`" + " ".join(cmd[1:])+"`...", MD.DEBUG)
        output = None
        try:
            output = subprocess.check_output(cmd)
        except subprocess.CalledProcessError as e:
            self.log("Process " + " ".join(cmd) + " returned non zero status " + str(e.returncode) + ". Abort.", MD.ERROR)
            if abort_on_failure:
                self.log(e.output, MD.FATAL)
            else:
                return e.returncode, e.output
        except:
            self.log("Unspecified failure on command "+" ".join(cmd) + ". Abort.", MD.FATAL)

        # 0 -- success code
        return 0, output

    def _run_cmd(self, cmd, abort_on_failure=True):
        return self.__run_cmd(map(str, cmd), abort_on_failure)

    def get_steps(self):
        step_files = sorted(glob.glob(self.get_run_dir() + "/run" + ("?"*len(self.__pattern % 0)) + ".dat"))
        steps = []
        prefix_len = len(self.get_run_dir() + "/run")
        suffix_len = len(".dat")
        for f in step_files:
            steps.append(int(f[prefix_len:-suffix_len]))
        return sorted(steps)

    def get_frames_in_step(self, step):
        if step < 0:
            step = self.get_current_step() + step+1
        suffix = self.__pattern % step
        bf = Binfile(self.get_run_dir() + "/" + "run" + suffix + ".dat")
        return len(bf.data)

    def save_state(self, filename=None):
        """
        Serialize MD state into marshal file
        :param filename:
        :type filename: str
        :return:
        """
        if filename is None:
            filename = MD.__dump_filename
        filename = os.path.abspath(filename)
        out = open(filename+".bak", 'wb')
        if out:
            self.log("Writing dump to %s "% filename)
            pickle.dump(self, out)
            out.close()
        else:
            self.log("Could not open file %s to write" % filename, MD.FATAL)
        shutil.move(filename+".bak",filename)

    def get_work_dir(self):
        return self.__work_dir;

    def __make_dir(self, dir_path):
        """ Creates directory if there is no file with the same name """
        self.log("Creating directory '"+dir_path+"'", MD.TRACE)
        if os.path.exists(dir_path):
            if os.path.isdir(dir_path):
                self.log("Directory '"+dir_path+"' already exists", MD.FATAL)
            else:
                self.log("Ordinary file '"+dir_path+"' already exists", MD.FATAL)
        else:
            os.makedirs(dir_path)

    def _set_hostname(self):
        self.__hostname = socket.gethostname()

    def __on_fatal(self):
        sys.exit(1)
