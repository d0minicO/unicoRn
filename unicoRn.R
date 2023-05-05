unicoRn <- function(base,
                    subs_name="unicoRn",
                    genes,
                    len="Full",
                    speciesToUse="xtropicalis|mdomestica|drerio|dmelanogaster|mmusculus|hsapiens|clfamiliaris|ggallus",
                    del_data=NULL,
                    check_delID=FALSE,
                    returnData=FALSE,
                    plotPdf=T){
  
  ## inputs
  # base = route / working directory
  # subs name = name of this list of genes to save in their own folder
  # genes = character string or vector eg "Gene", or c("Gene1", "Gene2")
  # len = the length of the first how many AAs to plot (or if you just want the whole length make len="whole" or any characters will do)
  # speciesToUse = character string of the species to use
  # del_data = location where the databse of IDs deleted from uniprot are
  # check_delID = TRUE or FALSE. Might be required if unable to pull some sequences from uniprot if they were deleted. Takes a long time so default should be false
  # returnData = TRUE or FALSE. If TRUE, then return a data.frame of the uniprot IDs, species, and sequences.
  
  #############
  # Libraries #
  #############
  
  require(msa)
  require(Biostrings)
  require(biomaRt)
  require(seqinr)
  require(tidyverse)
  require(data.table)
  require(tinytex)
  require(tools)
  
  
  #############
  # Functions #
  #############
  
  `%notin%` = function(x,y) !(x %in% y)
  
  makeChar = function(x) (c(x) %>% unlist %>% unname)
  
  
  ##########
  # Inputs #
  ##########
  
  dir.create(base, showWarnings = F)
  
  
  ###########
  # outputs #
  ###########
  
  alignDir = paste0(base, "alignments_", subs_name, "/")
  dir.create(alignDir, showWarnings = F)
  
  
  
  #################
  # Load the data #
  #################
  
  ### load the main ensembl mart
  ensembl <- useMart("ensembl")
  
  ## get the species specific datasets out
  ## first find which ones are there
  test = listDatasets(ensembl)
  
  ## load the species specific biomarts
  species_marts =
    test %>%
    filter(grepl(speciesToUse,dataset))%>%
    dplyr::select(dataset) %>%
    separate(dataset, into="ensembl_name",sep="_",extra = "drop",remove = F)
  
  
  # the list of all species that are supported
  ## the name ensembl uses and the name that a human wants to know
  all_species = data.frame(ensembl_name = c("xtropicalis",
                                            "mdomestica",
                                            "drerio",
                                            "dmelanogaster",
                                            "mmusculus",
                                            "hsapiens",
                                            "clfamiliaris",
                                            "ggallus"),
                           nice_name = c("Frog",
                                         "Opossum",
                                         "Zebrafish",
                                         "Drosophila",
                                         "Mouse",
                                         "Human",
                                         "Dog",
                                         "Chicken"))
  
  
  
  
  ## join these so that we get the correct order of human-readable species
  species = 
    left_join(species_marts,all_species,by="ensembl_name")$nice_name %>% 
    makeChar
  
  
  ########################
  # Get the AA sequences #
  ########################
  ## using biomart then uniprot
  
  
  ## loop through each mart
  data = data.frame()
  #species=c("Dog","Human")
  #spec=species[1]
  for(spec in species){
    
    
    # get the mart for this species
    mart = 
      left_join(species_marts,all_species,by="ensembl_name") %>% 
      filter(nice_name==!!spec) %>%
      dplyr::select(dataset) %>%
      makeChar
    
    
    message("Hello! loading all ensembl genes for ", spec, "\n")
    
    
    this_mart = useDataset(mart,mart=ensembl) # this is slow so only do once for the whole set of genes
    
    ### FETCH UNIPROT IDs USING BIOMART ###
    mapping <- getBM(
      attributes = c('uniprotswissprot', 'external_gene_name', "uniprotsptrembl"), 
      filters = 'external_gene_name',
      values = genes,
      mart = this_mart
    )
    
    # skip if nothing comes back from biomart
    if(nrow(mapping)==0){ 
      message("No records for given genes for ", spec, " moving on \n")
      next
    }
    ## loop through each gene and check it has been found
    ## only keep the "best" option
    ## in order of best to worst
    #1 uniprotswissprot >> Full curated uniprot entry
    #2 uniprotsptrembl nchar=6 >> Trembl entry (predicted)
    #3 uniprotsptrembl nchar=10 >> Trembl entry (predicted, obscure)
    
    df=data.frame()
    #gene=genes[1]
    for(gene in genes){
      
      
      message("Getting sequences for ", gene, "\n")
      
      # case insensitive match
      temp =
        mapping %>%
        filter(grepl(gene,external_gene_name, ignore.case = T))
      
      # have to replace any NAs with blank as otherwise this breaks the if statement below
      # make temp into data table to allow joining to the deleted IDs later
      temp =
        temp %>%
        mutate(across(everything(), ~ if_else(is.na(.), "", .)))
        data.table
      
      
      # if we didnt match that gene then skip it
      if (nrow(temp)==0){
        message("warning: no ", gene, " found in ", spec, " database \n")
        next
      }
      
      # if we only got one entry then that's cool
      if(nrow(temp)==1){
        message(gene, " was found uniquely ")
        
        # keep the uniprotID if there is one, or just keep the tremblID
        if(temp$uniprotswissprot!=""){
          message("taking uniprot ID \n")
          id = makeChar(temp$uniprotswissprot)
          
        }
        
        else if(nchar(temp$uniprotsptrembl)==6){
          message("taking TREMBL ID short \n")
          id = makeChar(temp$uniprotsptrembl)
          
        } else if (nchar(temp$uniprotsptrembl)==10){
          message("taking TREMBL ID long \n")
          id = makeChar(temp$uniprotsptrembl)
          
        } else {
          warning(gene, " encountered a problem, check it \n")
        }
        
      } else if (nrow(temp)>1){ # end of the unique returns
        message(gene, " has more than one match ")
        
        
        # look for the uniprot ID
        temp1 =
          temp %>%
          filter(uniprotswissprot!="") %>%
          dplyr::slice(1)
        
        # look for the short TREMBL ID
        temp2 =
          temp %>%
          filter(nchar(uniprotsptrembl)==6) %>%
          dplyr::slice(1)
        
        # look for the long TREMBL ID
        temp3 =
          temp %>%
          filter(nchar(uniprotsptrembl)==10) %>%
          dplyr::slice(1)
        
        
        
        if(nrow(temp1)==1){
          message("taking uniprot ID \n")
          id = makeChar(temp1$uniprotswissprot)
          
        } else if (nrow(temp2)==1){
          message("taking short TREMBL ID \n")
          
          id =
            temp2 %>%
            filter(nchar(uniprotsptrembl)==6)%>%
            dplyr::select(uniprotsptrembl) %>%
            makeChar
          
        } else if (nrow(temp3)==1){
          message("taking long TREMBL ID \n")
          id =
            temp3 %>%
            filter(nchar(uniprotsptrembl)==10)%>%
            dplyr::select(uniprotsptrembl) %>%
            makeChar
          
        } else {
          message(" something went WRONG \n")
        }
        
        
      } # end of the multi row returns
      
      
      
      temp_df = data.frame(Species = spec,
                           Gene = gene,
                           Uniprot_id = id)
      
      df=rbind.data.frame(df,temp_df)
      
    } # end of for gene loop
    
    
    # save into the main df
    data = rbind.data.frame(data,df)
    
  }
  
  
  ## should we check the deleted IDs database?
  if(check_delID){
    ## now once all the uniprot IDs are in we can remove any deleted IDs
    ### all accessions deleted by uniprot
    ## huge table takes ages to load
    ## only work on in data.table syntax not dplyr!!
    message("Now loading the uniprot deleted IDs database which will take a minute or two...\n")
    del_accs = readRDS(del_data)
    # set keys to help join faster
    setkey(del_accs, "ID")
    
    # make the ids table into a data table to allow quick join to the deleted IDs
    data =
      data %>% data.table()
    # set the key
    setkey(data,"Uniprot_id")
    
    deleted = makeChar(data[del_accs, nomatch = 0])
    
    # remove any deleted IDs
    data =
      data %>%
      filter(Uniprot_id %notin% deleted)
    
    if(length(deleted)==0){
      message("No deleted IDs to remove! \n")
    } 
    else if (length(deleted)!=0) {
      message("Successfully removed ", length(deleted), " deleted IDs \n")
    } 
  } else if (!check_delID){
    message("Not checking uniprot for deleted IDs to save time, skipping \n")
  }
  
  # now loop through each row of the output table and extract the sequence from uniprot
  message("Now pulling ", nrow(data), " sequences from uniprot \n")
  df = data.frame()
  for(i in 1:nrow(data)){
    
    temp_seq = "Empty"
    
    spec = data[i,1]
    gene = data[i,2]
    id = data[i,3]
    
    
    ### get the AA sequence from uniprot
    acc_url = paste0("https://www.uniprot.org/uniprot/",id,".fasta")
    temp_seq=paste0(read.csv(url(acc_url))[,1],collapse="")
    
    ## set up a useful warning message
    if(temp_seq == "Empty"){
      warning("Uniprot sequence not found. Try running with check_delID=TRUE and provide location of del_data \n")
      break
    }
    
    
    # piece together the df
    temp_df = data.frame(gene,id,temp_seq,spec)
    df = rbind.data.frame(df,temp_df)
  } # end of the for each id loop
  
  
  message("Done gathering sequences from uniprot")
  
  
  # set up the levels of the species in the df correctly
  ## this makes sure the alignment is done from Human down to Fly
  data = as_tibble(df)
  
  data =
    data %>%
    dplyr::rename(Species=spec,
                  Gene=gene)
  
  data$Species = factor(data$Species, levels=c("Human",
                                               "Mouse",
                                               "Dog",
                                               "Opossum",
                                               "Chicken",
                                               "Frog",
                                               "Zebrafish",
                                               "Drosophila"))
  
  
  
  data =
    data %>% arrange(Species)
  
  ## get a new column that has a concise conjoined name
  data = 
    data %>% 
    mutate(name=paste0(Gene,"_",Species))
  
  
  ## should we return the data?
  if(returnData){
    message("Just returning the data for you and not aligning or plotting ")
    return(data)
  } else if (!returnData){
    
    message("Moving on to alignment \n")
    
    ##############
    # ALIGNMENTS #
    ##############
    
    #g=genes[3]
    tex_files = data.frame()
    for (g in genes){
      
      # outputs
      alFile = paste0(alignDir,"UnicoRn_analysis_",g,".fasta")
      texFile = paste0(alignDir,"UnicoRn_analysis_",g,".tex")
      pdfFile = paste0(alignDir,"UnicoRn_analysis_",g,".pdf")
      
      # filter on just this gene
      temp =
        data %>%
        filter(Gene==g)
      
      
      if(nrow(temp)<2){
        message("Not enough sequences for ",g, " to align, skipping")
        next
      }
      
      # get the AAStringSet object and give it proper names
      seqs = AAStringSet(temp$temp_seq)
      names(seqs) = temp$name
      
      
      ## subset the first X (len) AAs
      if(is.numeric(len)){
        seqs = subseq(seqs,1,len)
      }
      
      # do the alignment
      aligned = msa(seqs,order = "input")
      
      # save the alignment in a way texshade can compile
      msaPrettyPrint(aligned,
                     alFile = alFile,
                     #file = texFile,
                     output="asis",
                     #showNames="none",
                     #showLogo="top",
                     #logoColors="rasmol",
                     #shadingMode="similar",
                     #showLegend=FALSE,
                     askForOverwrite=FALSE#,
                     #furtherCode=c("\\defconsensus{.}{lower}{upper}",
                     #             "\\showruler{1}{top}")
      )
      
      
      message("Saved fasta alignments for ", g)
      
      
      ## save the tex file
      ### modify the display options as you please!
      #### best if you know some latex :)
      writeLines(c(
        "\\documentclass[preview]{standalone}",
        "\\usepackage{texshade}",
        "\\usepackage{inconsolata}",
        "\\usepackage{geometry}%[showframe]",
        
        "\\begin{document}",
        paste0("\\begin{texshade}{",alFile,"}"),
        # IDENTITY HIGHLIGHTING
        "\\shadingmode[allmatchspecial]{identical}",
        "\\nomatchresidues{Gray70}{White}{upper}{bf}",
        "\\conservedresidues{Black}{LightCyan}{upper}{bf}",
        "\\allmatchresidues{White}{Red}{upper}{bf}",
        
        # CHEMICAL HIGHLIGHTING
        #\\shadingmode[chemical]{functional}
        
        # STRUCTURAL HIGHLIGHTING
        #\\shadingmode[structure]{functional} 
        
        # HIGHLIGHTING THRESHOLD
        "\\threshold[100]{50}", #% [high]{low}
        
        # CAPTION - textbf{} specifies bold text
        #\\showcaption[bottom]{\\textbf{Protein MSA with Similarity Highlighting}} 
        
        # LEGEND
        "\\showlegend",
        "\\movelegend{0cm}{0cm}", #% {Horizontal}{Vertical} offsets
        
        # TOP NUMBERING
        "\\showruler{1}{top}",
        #\\hidenumbering
        
        # SIDE NUMBERING
        #"\\shownumbering",
        
        # CONSENSUS
        #\\hideconsensus
        "\\showconsensus[ColdHot]{bottom}",
        "\\defconsensus{.}{lower}{upper}",
        
        # FINGERPRINTING
        #\\fingerprint{200}
        
        # SEQUENCE LOGO
        #"\\showsequencelogo{top}",
        #\\showlogoscale{leftright}
        #\\dofrequencycorrection
        #\\setends{1}{1..200}
        
        # TEXT SIZE (see README)
        "\\namesfootnotesize",
        "\\residuesfootnotesize",
        "\\legendfootnotesize",
        "\\numberingtiny",
        
        # end document
        "\\end{texshade}",
        "\\end{document}"),
        texFile)
      
      message("Printing the pdf for you now... \n")
      
      
      ## dont plot PDF if not requested
      ## can help debug buggy texi2pdf
      if(plotPdf){
        ## save pdf in root folder
        tools::texi2pdf(texFile, clean=TRUE)
        #rename the temp file to the desired output
        pdf_file_temp = paste0("~/","UnicoRn_analysis_",g,".pdf")
        file.rename(pdf_file_temp, pdfFile)
        
      }

    }
    
  }
  
  
  
} # end of function
