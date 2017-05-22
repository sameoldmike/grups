#                           $1                  $2              $3     
# usage:    thisscript.sh   input_filename_base input_data_dir  output_dir_basename
# ex:       thisscript.sh LeipzigSims_c0.01,0.01,0.01,0.01,0.01,0.01_cov5.0 ./pedigree_sims_complete.2015-02-05.py.LeipzigSims/ ConcatenatedReps/

cat $2/*$1*.out            >  $2$3$1.out
cat $2/*$1*.numvars        >  $2$3$1.numvars
cat $2/*$1*.numSNPscov     >  $2$3$1.numSNPscov
cat $2/*$1*.pwexp          >  $2$3$1.pwexp                                    
cat $2/*$1*.pwexp_sibs     >  $2$3$1.pwexp_sibs    
cat $2/*$1*.pwexp_gpgc     >  $2$3$1.pwexp_gpgc    
cat $2/*$1*.pwexp_twins    >  $2$3$1.pwexp_twins    
cat $2/*$1*.pwexp_fcous    >  $2$3$1.pwexp_fcous
cat $2/*$1*.labels         >  $2$3$1.labels    
