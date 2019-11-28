

def get_trjtool_input(coord_dict):
    template_input = """############################################################
#                  DCD Trajectories                        #
############################################################
#
#------ filename of topology file for trajectories
#------ currently only supports pdb and charmm psf format
fnpsf  {reference_structure}

#------ filename of trj files, start & end file index
#------ currently only supports dcd and netcdf format
fntrj  {netcdf_pattern}      {netcdf_number} {netcdf_number}

#------ <stride>  <first>  <last>
#------ stride and first&last frame# for each trj file
#------ if <first> <first> are absent: 1 ~ LEN
#stride 10  1  5000
stride  {frame_stride} {frame_first_number} {frame_last_number}



############################################################
#                    Frame Alignment                       #
############################################################
#
#====== whether to superimpose frames
align  false  #<=== true/false

#------ which segments to be focused. Using $ for blank
#------ empty means all segments
#seg_id  SEGA SEGB
seg_id  $$$$

#------ which resids to search; blank means all resids
resid  7-11 14 18-60 66 68 71-73

#------ which residue names to search; blank means all names
resname

#------ atoms to be fitted. Normally 'N CA C' or 'CA'
atname  CA

#------ which file name to save the rmsd data
fnout  rmsd.tab

#------ output as ascii(true) or binary(false)
ascii  true



############################################################
#                     Job Definition                       #
############################################################
#
#
#>>>>>>>>>>>>> Calculating radius of gyration <<<<<<<<<<<<<#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rgyr  false  #<=== true/false

#------ which segments to be analyzed. Using $ for blank
#------ empty means all segments
#seg_id  PROT PEP$ XWAT
seg_id  $$$$

#------ if leave it blank, all resid are analyzed
#resid  1-5 7-21 26 30-76
resid

#------ if leave it blank, all resname are analyzed
#resname  LEU ILE VAL
resname

#------ use atom # if you want to fix this atom
#------ e.g. in the case of spin label
#atname  HN  949
atname  heavy

#------ which file name to save the result
fnout  rgyr.tab

#------ output as ascii(true) or binary(false)
ascii  true
#===========================================================


#>>>>>>>>>>>>>>>   Extracting coordinates   <<<<<<<<<<<<<<<#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
coor  true  #<=== true/false

#------ which segments to be analyzed. Using $ for blank
#------ empty means all segments
#seg_id  PROT PEP$ XWAT
seg_id  {segment_ids}

#------ if leave it blank, all resid are analyzed
#resid  1-5 7-21 26 30-76
resid  {residue_ids}

#------ if leave it blank, all resname are analyzed
#resname  LEU ILE VAL
resname {residue_names}

#------ use atom # if you want to fix this atom
#------ e.g. in the case of spin label
#atname    HN  949
atname  {atoms}

#------ which file name to save the result
fnout  {output_filename}

#------ output as ascii(true) or binary(false)
ascii  {is_ascii}
#===========================================================


#>>>>>>>>>>>>>>>     Extracting vectors     <<<<<<<<<<<<<<<#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
vect  false  #<=== true/false

#------ which segments to be analyzed. Using $ for blank
#------ empty means all segments
#seg_id  PROT PEP$ XWAT
seg_id  $$$$

#------ if leave it blank, all resid are analyzed
#resid  1-5 7-21 26 30-76
resid

#------ if leave it blank, all resname are analyzed
#resname  LEU ILE VAL
resname

#------ use atom # if you want to fix this atom
#------ e.g. in the case of spin label
#atname  HN  949
atname  N HN

#------ which file name to save the result
fnout  vect.dat

#------ output as ascii(true) or binary(false)
ascii  false
#===========================================================


#>>>>>>>>>>>>>>>    Extracting distances    <<<<<<<<<<<<<<<#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dist  false  #<=== true/false

#------ which segments to be analyzed. Using $ for blank
#------ empty means all segments
#seg_id  PROT PEP$ XWAT
seg_id  $$$$

#------ if leave it blank, all resid are analyzed
#resid  1-5 7-21 26 30-76
resid

#------ if leave it blank, all resname are analyzed
#resname  LEU ILE VAL
resname

#------ use atom # if you want to fix this atom
#------ e.g. in the case of spin label
#atname    HN  949
atname  N HN

#------ which file name to save the result
fnout  dist.tab

#------ output as ascii(true) or binary(false)
ascii  true
#===========================================================


#>>>>>>>>>>>>>>>     Extracting angles      <<<<<<<<<<<<<<<#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
angl  false  #<=== true/false

#------ which segments to be analyzed. Using $ for blank
#------ empty means all segments
#seg_id  PROT PEP$ XWAT
seg_id  $$$$

#------ if leave it blank, all resid are analyzed
#resid  1-5 7-21 26 30-76
resid

#------ if leave it blank, all resname are analyzed
#resname  LEU ILE VAL
resname

#------ use atom # if you want to fix this atom
#------ e.g. in the case of spin label
#atname    HN  949
atname  HN N CA

#------ which file name to save the result
fnout  angl.tab

#------ output as ascii(true) or binary(false)
ascii  true
#===========================================================


#>>>>>>>>>>>>>>> Extracting dihedral angles <<<<<<<<<<<<<<<#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dihe  false  #<=== true/false

#------ which segments to be analyzed. Using $ for blank
#------ empty means all segments
#seg_id  PROT PEP$ XWAT
seg_id  $$$$

#------ if leave it blank, all resid are analyzed
#resid  1-5 7-21 26 30-76
resid

#------ if leave it blank, all resname are analyzed
#resname  LEU ILE VAL
resname

#------ use atom # if you want to fix this atom
#------ e.g. in the case of spin label
#atname  HN  949
#------ use +/- to indicate next/previous residue
#------ e.g. to get backbone phi and psi, use
#atname  -C N CA C    N CA C +N
atname  HN N CA HA

#------ which file name to save the result
fnout  dihe.tab

#------ output as ascii(true) or binary(false)
ascii  true
#===========================================================
"""
    return template_input.format(**coord_dict)
