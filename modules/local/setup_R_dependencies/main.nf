/*
   This installs a couple R packages that are not included in the
   base tidyverse singularity image we are using.
*/
process SETUP_R_DEPENDENCIES {
  label      'process_low'
  tag        "${R_lib_dir}"
  publishDir "${params.outdir}"


  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.2.2"
  }

  input:
  path (install_script)   // the path to a script that will handle package installation
  path (R_tar_dir)        // the path to tar.gz files of R packages to be installedj
  // val  (R_lib_dir)        // the name of a dir that will be created to contain newly installed
                          // R pkgs.  This not a path type because not an existing path - will be created

  output:
  val  (true)           , emit: setup_complete // a flag to indicate this is complete
  path ("R_lib_dir")    , emit: R_lib_dir
  path "versions.yml"   , emit: versions

  script:

  if (workflow.containerEngine == 'singularity') {
  """
    mkdir -p R_lib_dir
    Rscript ${install_script} ${R_tar_dir} R_lib_dir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
  """
  } else {
  """
    echo "setup not necessary for conda environment"
    # just create an empty directory 
    mkdir R_lib_dir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
  """
  }
}

