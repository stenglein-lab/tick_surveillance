include { BLAST_MAKEBLASTDB } from '../../../modules/nf-core/blast/makeblastdb/main'                                                                                                                                                       

// MAKE_BLAST_DB options
process {
    withName: 'BLAST_MAKEBLASTDB' {

        ext.args   = [
                "-dbtype nucl"
            ].join(' ').trim()
        ]
    }
}

/*
   Setup indexes and dictionaries needed by downstream processes.
*/
workflow SETUP_INDEXES {

  take: 
  refseq_fasta

  main:
  BLAST_MAKEBLASTDB(refseq_fasta)
  
  emit:
  db       = BLAST_MAKEBLASTDB.out.db
  versions = BLAST_MAKEBLASTDB.out.versions
}

