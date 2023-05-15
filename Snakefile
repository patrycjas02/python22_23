assemblies, = glob_wildcards('/reserve/data/B2/genomes/{assembly}.fna')
    
rule all:
   input:
       'output/tree.txt'
        
rule blastn:
   params:
       outfmt = 6,
       cols   = 'sseqid sseq qcovs evalue'
   input:
       query   = 'input/rpoB.fna',
       subject = '/reserve/data/B2/genomes/{assembly}.fna'
   output:
       'output/blastn/blastn_{assembly}.tsv'
   log:
       'log/blastn/blastn_{assembly}.log'
   shell:
       '''blastn -outfmt     "{params.outfmt} {params.cols}" \
                 -query    {input.query}                          \
                 -subject  {input.subject}                          \
                 -out      {output}                          \
                  > {log} 2>&1
       '''
        
        
rule filter:
    params:
        qcovs = 100,
        cols  = rules.blastn.params.cols
    input:
        rules.blastn.output
    output:
        'output/filter/filter_{assembly}.tsv'
    log:
        'log/filter/filter_{assembly}.log'
    shell:
        '''scripts/filter.py --qcovs    {params.qcovs} \
                             --cols    "{params.cols}" \
                             --input    {input}        \
                             --output   {output}       \
                               > {log} 2>&1
        '''
    
    
rule convert:
    input:
        rules.filter.output
    output:
        'output/convert/convert_{assembly}.fna'
    log:
        'log/convert/convert_{assembly}.log'
    shell:
        '''scripts/convert.py --input	{input}			\
							  --output	{output}		\
        '''
    
    
rule merge:
    params:
        mask = rules.convert.output[0].replace('{assembly}', '*')
    input:
        expand(rules.convert.output, assembly=assemblies)
    output:
        'output/merged.fna'
    log:
        'log/merge.log'
    shell:
        'cat {params.mask} > {output} 2> {log}'
        
        
rule clustalo:
    params:
        infmt  = 'fasta',
        outfmt = 'fasta'
    threads:
        8
    input:
        rules.merge.output
    output:
        'output/clustalo.fna'
    log:
        'log/clustalo.log'
    shell:
        '''clustalo --infmt    {params.infmt} \
                    --outfmt   {params.outfmt} \
                    --threads  {threads} \
                    --infile   {input} \
                    --outfile  {output} \
                      > {log} 2>&1
        '''
        
tree_name = 'rpoB'
        
rule raxml:
    threads:
        8
    input:
        rules.clustalo.output
    output:
        'output/RAxML_bestTree.' + tree_name
    log:
        'log/raxml.log'
    params:
        model  = 'GTRCAT',
        outdir = os.getcwd() + '/output',
        name   = tree_name
    shell:
        '''raxml -p 1               \
                 -T {threads}       \
                 -s {input}         \
                 -w {params.outdir} \
                 -n {params.name}   \
                 -m {params.model}  \
                  > {log} 2>&1
        '''
        
rule drawtree:
    input:
        rules.raxml.output
    output:
        'output/tree.txt'
    log:
        'log/drawtree.log'
    shell:
        '''scripts/drawtree.py --input  {input}  \
                               --output {output} \
                                 > {log} 2>&1
        '''

