#test



stepwise_valued_simulation <- function(name_of_outputfile = "random",
                                 save_directory = "current",
                                 durchlaeufe_ = 10000, 
                                 mrna_ = "mRNA.fa",
                                 srna_ = "sRNA.fa",
                                 translation_start_ = 34,
                                 core_ = 3,
                                 time_steps = 5,
                                 wiederholungen = 1,
                                 target_grenze_ = -18,
                                 store_ = 50,
                                 target_value_threshold_ = 0.8,
                                 srna_value_threshold_ = 0.8,
                                 target_constrain_ = TRUE,
                                 rna_constrain_ = TRUE,
                                 accessibility = TRUE,
                                 keep_all_rnas = FALSE,
                                 change_start_values = 0.5,
                                 percentage_of_pos_values = 0.5,
                                 start_target_number = 1){
  ##### Beschreibung #####
  
  #target_value_threshold_: Target Value Constrain. 0.5 = 50%. Makes sure that the old Target Value wont decrease to much after each mutation checkpoint.
  #srna_value_threshold_ : Complete Value Constrain. 0.5 = 50%. Function compares Old sRNA value to new sRNA value and decide based on this parameter if the Mutations will be allowed.
  # target_grenze_: Target threshold. Defines a Energie threshold. All mRNAs wich bounding Energy to the sRNA is below this threshold, will be seen as Targets.
  #change_start_values: Defines the Value of the First sRNA Target based on the mean of all positive mRNA Values. 1.2 = 120% of the mean (WARNING: 10 = max value not 1000%)
  # percentage_of_pos_val: Percentage of mRNAs with a positive value. 0.4 = 40% of the mRNAs will have a positive Value
  # start_target_number_: Warning! this one can be irritating. Defines the number of Start Targets which will have a positive Value. All other Targets will have the Value 0.
  ## If the sRNA only has 4 Targets in the beginning, then 4 is the max you can type in here.
  #mrna_name: Data name of the mRNA file. Must be a RData file with the listed mRNA`s and their values as attributes. (See function value_for_mrna)
  #durchlaeufe: Numbers of Mutations per mutation cycle
  #srna: sRNA file (.fa or .txt but fasta format)
  #core: specify the numbers of Cores that will be used. (Note: More Cores -> faster calculation but sometimes there occur errors if to many cores are used and you get 0byte Outputs)
  #wiederholungen: times the Function will repeat itself. core * wiederholungen = number of trials
  #time_steps: number of Mutation cycles that will be preformed. time_steps * durchlaeufe = Complete number of Mutations that will be inserted
  #store: if 200 the function will check the sRNa and mRNA`s after 200 inserted Mutation and adopt or reject them, based on the Target Values.
  #target_constrain / rna_constrain: FALSE equals target_con / value_con = 0. 
  # accessibility: TRUE = intaRNA call will include the accessibility of the RNAs
  # keep_all_rnas: FALSE = RNAs wont be safed in the output.
  ##### Beschreibung Ende #####
  
  
  ##### Packages #####
  require(seqinr)
  require(doMC)
  
  ##### Packages Ende #####
  if (target_value_threshold_ > 1 | srna_value_threshold_ > 1){
    stop("threshold  too big")
  }
  
  
  ##### Funktionen #####
  srna_value_evolution <- function(srna,
                                   mrna,
                                   translation_start,
                                   durchlaeufe,
                                   target_grenze = -13,
                                   store = 50,
                                   target_value_threshold = 0.8,
                                   srna_value_threshold = 0.8,
                                   target_constrain = TRUE,
                                   rna_constrain = TRUE){
    #srna/mrna: sRNA and mRNAs in lists 
    #target_value_threshold: Target Value Constrain. 0.5 = 50%. Makes sure that the old Target Value wont decrease to much after each mutation checkpoint.
    #srna_value_threshold : Complete Value Constrain. 0.5 = 50%. Function compares Old sRNA value to new sRNA value and decide based on this parameter if the Mutations will be allowed.
    # target_grenze_: Target threshold. Defines a Energie threshold. All mRNAs wich bounding Energy to the sRNA is below this threshold, will be seen as Targets.
    #durchlaeufe: Numbers of Mutations per mutation cycle
    #store: if 200 the function will check the sRNa and mRNA`s after 200 inserted Mutation and adopt or reject them, based on the Target Values.
    #target_constrain / rna_constrain: FALSE equals target_con / value_con = 0. 

    
    #BESCHREIBUNG OUTLIST!!
    # Die Outlist ist ein Listen ELement mit 7 eintraegen:
    # 1.) die aktuelle sRNA
    # 2.) die aktuellen mRNA's
    # 3.) anzahl an Mutationen die tatsächlich in die sRNA/mRNAs einefügt wurden
    # 4.) Die Positionsinformationen der aktuellen targets
    # 5.) Die ursprüngliche sNRA
    # 6.) Die Ursprünglichen Targets


    
    ######## Packages ########
    require(seqinr)
    ######## Packages Ende ########
    
    ######## Functions ########
    stack_mutations <- function(mutant_list1, translation_start1, store, the_rna_constrain){
      #Input: mutant_list1: Eine kopie der Outlist in dem sich alle Informationen der RNA's befinden. Wichtig sind hier eigentlich nur die ersten 2 Eintraege
      #                     1.) die sRNA 2.) alle mRNAs
      #       translation_start1: Position des ersten Nucleotides nach dem start Codon
      #       store: anzahl der Mutationen die Gesammelt werden sollen.
      #       the_rna_constrain: TRUE oder FALSE logical der Angibt ob die mRNA's in ihrer Mutation restriktiert sind. Wenn FALSE, können alle mrNAs Komplett zufaellig mutieren 
      #     
      #
      #Diese Funktion fügt eine beliebige Anzahl an Mutationen in die bestehenden RNAs ein. Nach jeder Mutation passieren folgende dinge:
      #   1.) Check ob die Mutation die Sequenz veraender hat, falls nicht wird die Mutation verworfen
      #   2.) Falls die Mutation nicht verworfen wurde, wird sie in einer Liste gespeichert, und der Name der listenposition enspricht der Position der RNA in der Gesamt RNA Liste
      
      #       Falls eine RNA 2 mal Mutiert, ersetzt sie die alte RNA in der Liste -> so sind mehrere Mutation in einer RNA möglich
      #       Es wird sicher gestellt, das sich eine Mutierte sRNA IMMER am anfang der Mutations Liste befindet, um sie so spaeter von den anderen Mutationen Trennen zu können
      
      # Output: Eine Liste mit 2 Eintraegen: 1.) Die Liste mit den Mutierten RNA's (wichtig für die mRNA Kontrolle)
      #                                     2.) Die mutant_list1 in denen sich alle sequenzen befinden (auch die Mutierten) (Wichtig für die sRNA Kontrolle)
      
      temp_mutation_list1 <- list()
      count1 <- 0
      for(i in 1:store){
        current_mutation_info <- onemutation(mutant_list1[[1]], mutant_list1[[2]], translation_start1, the_rna_constrain)
        
        #Kontrolle ob die Mutation eine Veraenderung hervorgerufen hat
        if(current_mutation_info[[3]] == 0){ 
          next
        }
        count1 <- count1 +1
        #Kontrolle Ob die Mutation in einer mRNA stattfand (falls die sRNA mutiert ist, dann ist current_mutation_info[[2]] = 0)
        if(current_mutation_info[[2]] != 0){
          #Mutierte RNA ersetzt alte RNA in der mutant_list
          mutant_list1[[2]][[current_mutation_info[[2]]]] <- current_mutation_info[[1]]
          
          #RNA wird der temp-muation_list hinzugefügt
          if(length(which(names(temp_mutation_list1) == paste(current_mutation_info[[2]], collapse = ""))) == 0){
            temp_mutation_list1[[length(temp_mutation_list1) + 1]] <- current_mutation_info[[1]]
            names(temp_mutation_list1)[length(temp_mutation_list1)] <- paste(current_mutation_info[[2]], collapse = "")
          } else {
            temp_mutation_list1[[paste(current_mutation_info[[2]], collapse = "")]] <- current_mutation_info[[1]]
          }
          
        } else{
          #Mutation fand in der sRNA statt
          mutant_list1[[1]] <- current_mutation_info[1]
          if(length(which(names(temp_mutation_list1) == paste(current_mutation_info[[2]], collapse = ""))) == 0){
            temp_mutation_list1 <- c(list(current_mutation_info[1]), temp_mutation_list1)
            names(temp_mutation_list1)[1] <- current_mutation_info[[2]]
          } else {
            temp_mutation_list1[[1]] <- current_mutation_info[1]
          }
        }
      }
      return(list(mutant_list1, temp_mutation_list1, count1))
    }
    
    onemutation <- function(srna, trna, first_coding_nuc, the_rna_constrain) {
      ## Input: srna: eine sRNA 
      #         trna: eine liste mit mrna's
      #         first_coding_nuc: Position des ersten Nucleotides nach dem Start-Codon
      #         the_rna_constrain: TRUE oder FALSE logical der Angibt ob die mRNA's in ihrer Mutation restriktiert sind. Wenn FALSE, können alle mrNAs Komplett zufaellig mutieren (ohne unten genannte Einschraenkungen).
      #
      
      # Diese Funktion fügt eine zufaellige Mutation in eine Zufaellige rna ein, dabei ist die Mutationswahrscheinlichkeit für jede RNA abhaengig von ihrer laenge(sRNA hat eine niedrigere Mutationswahrscheinlichkeit wenn sie kürzer ist)
      # Beim einfügen einer Mutation werden folgende Dinge berücksichtigt:
      
      #Findet die Mutation in der mRNA statt, dann werden die Mutationswarscheinlichkeiten der UTR's aus einer Pssm-Matrix gelese 
      #DIe PSSM-MATRIX MUSS IM GLOBAL VORHANDEN SEIN / Name: utr_pssm (siehe Funktion make_pssm_utr)
      
      #Ist die Mutation im Coding bereich, dann wird ermittelt um welches Triplet bzw Aminosaeure es sich handelt. Danach
      #wird anhand einer wahrscheinlichkeit für das einsetzen einer neuen AS mutiert. Diese werden aus einer PAM Matrix abgelesen
      #ACHTUNG PAM MATRIX MUSS IM GLOBAL VORHANDEN SEIN Name: pam (siehe Funktion initial_pam)
      #ACHTUNg es muss auch die Matrix ttoas im Global vorhanden sein, diese Übersezt Triplets in die Zugehörigen AS (siehe Funktion initial_triplet_to_as)
      
      #ACHTUNG: falls sich die Mutation am ende der Sequenz befindet und es z.b einen Überhang von 2 Nucleotiden gibt (Es sich also aufgrund eines Fehlenden Nucleotides kein TRiplet bilden laesst)
      # Dann wird die Mutation in so einem Überhang NICHT erlaubt
      #Die sRNA kann frei mutieren und hat keine einschraenkungen
      #
      ##Output: Eine Liste mit 3 Eintraegen: 1.) Mutierte RNA
      #                                    2.) Position der Mutierten RNA in der eingegebenen RNA liste (falls 0, ist die sRNa Mutiert)
      #                                    3.) Information ob die Mutation die Nucleotidsequenz veraendert hat (entweder 0 oder 1)
      #                                        0: Die Sequenz ist gleich geblieben / 1: Die sequenz hat sich veraendert
      
      
      #Initialisiert ein paar variablen der Funktion EIN mal! Achtung Variablen sind nun im Global vorhanden 
      if(length(which(ls() == "first_onemutation"))==0){
        len_trna_onemutation <<- length(trna[[1]])
        len_srna_onemutation <<- length(srna[[1]])
        len_all_onemutation <<- length(srna) + length(trna) *length(trna[[1]])
        first_onemutation <<- "Generated"
      }
      
      mutant_info <- vector("list", 3)
      names(mutant_info) <- c("mutant","rna pos(if 0 -> srna)","is new?")
      #zu Mutierende RNA wird ermittelt
      choose <- sample(1:len_all_onemutation, 1)
      rna <- ceiling((choose - len_srna_onemutation) / len_trna_onemutation)
      if(rna < 0){
        rna <- 0
      }
      
      if(rna != 0){
        selectet_rna <- trna[[rna]]
        pos <- sample(1:len_trna_onemutation, 1)
      } else {
        selectet_rna <- srna
        pos <- sample(1:len_srna_onemutation, 1)
      }
      #rna: 0 wenn es sRNA getroffen hat oder Position der mRNA
      #pos: Nucleotidposition
      mutant_info[[2]] <- rna
      if((the_rna_constrain == FALSE) & (rna != 0)){
        nucleotide <- c("a","g","c","t")
        mut_nuc <- sample(nucleotide,1)
        mutant <- selectet_rna
        mutant[[pos]] <- mut_nuc
        mutant_info[[1]] <- mutant
        if(paste(mutant, collapse = "") == paste(selectet_rna, collapse = "")){
          mutant_info[[3]]<- 0
        } else {
          mutant_info[[3]]<- 1
        }
        return(mutant_info)
      }
      
      if(rna == 0){
        nucleotide <- c("a","g","c","t")
        mut_nuc <- sample(nucleotide,1)
        mutant <- srna[[1]] #Verursacht vlt Probleme
        mutant[[pos]] <- mut_nuc
        mutant_info[[1]] <- mutant
        if(paste(mutant, collapse = "") == paste(srna[[1]], collapse = "")){
          mutant_info[[3]]<- 0
        } else {
          mutant_info[[3]]<- 1
        }
      } else if(pos >= first_coding_nuc){
        mutant_info_zip <- coding_reg_mutation(selectet_rna, first_coding_nuc, pos)
        mutant_info[[1]] <- mutant_info_zip[[1]]
        mutant_info[[3]] <- mutant_info_zip[[2]]
      } else{
        mut_nuc <- sample(rownames(utr_pssm),1,prob = utr_pssm[,pos])
        mutant <- selectet_rna
        mutant[pos]<-mut_nuc
        mutant_info[[1]] <- mutant
        if (paste(mutant, collapse = "") == paste(selectet_rna, collapse = "")){
          mutant_info[[3]] <- 0
        } else{
          mutant_info[[3]] <- 1
        }
      } 
      #mutant_info besteht aus einer liste mit der Mutierten Sequenz und informationen darüber ob die sequez die srna ist, sowie ob die Mutation was veraendert hat.  
      
      return(mutant_info)
      
    }
    
    coding_reg_mutation <- function(sequence_complete, first_coding_nuc, mut_position){
      #Input: sequence_complete: eine Nucleotid Sequenz die zu Mutieren ist
      #       first_coding_nuc: Position des ersten Nucleotides nach dem start Codon
      #       mut_position: Position an dem die Mutation in de rSequenz eingefügt werden soll
      
      #AcHTUNG benötigt ttoas matrix und PAM- MAtrix
      #ttoas: Matrix welche die Triplets in AS übersetzt
      #Extrahiert aus der Position das Triplet, welches Mutiert wird und Mutiert das Tripletnach Wahrscheinlichkeiten aus der PAM Matrix
      #ACHTUNG: falls sich die Mutation am ende der Sequenz befindet und es z.b einen Überhang von 2 Nucleotiden gibt (Es sich also aufgrund eines Fehlenden Nucleotides kein TRiplet bilden laesst)
      # Dann wird die Mutation in so einem Überhang NICHT erlaubt
      #
      #Output: eine Liste mit 2 Eintraegen: 1.) Die Mutierte Sequenz
      #                                    2.) Information ob die Mutation die Nucleotidsequenz veraendert hat (entweder 0 oder 1)
      #                                        0: Die Sequenz ist gleich geblieben / 1: Die sequenz hat sich veraendert 
      
      if (first_coding_nuc > mut_position){
        print("Error: Mutation liegt nicht in der coding Region!  Funktion: coding_reg_mutation")
        stop("Error: Mutation liegt nicht in der coding Region!  Funktion: coding_reg_mutation")
      } else{
        mut_position2 <- mut_position - first_coding_nuc + 1
        len_coding_reg <- length(sequence_complete) + 1 - first_coding_nuc
        ueberstand <- len_coding_reg %% 3
        if( mut_position > length(sequence_complete)){
          print("ERROR: Mutations Position liegt ausserhalb der Sequenz Funktion: coding_reg_mutation ")
          stop("ERROR: Mutations Position liegt ausserhalb der Sequenz Funktion: coding_reg_mutation ")
        }
        if((len_coding_reg - mut_position2) >= ueberstand){
          pos_triplet <- ceiling(mut_position2/3)
          pos_nucleotid <- mut_position2%%3
          if(pos_nucleotid == 0){
            pos_nucleotid <- 3
          }
          mut_triplet <- sequence_complete[(((pos_triplet - 1) * 3) + first_coding_nuc):((pos_triplet * 3) -1 + first_coding_nuc)]
          possible_triplets <- vector("list",4)
          nucleotide <- c("a","g","c","t")
          for(i in 1:4){
            possible_triplets[[i]]<-replace(mut_triplet,pos_nucleotid,nucleotide[i])
          }
          po_tr <-c()
          for(i in 1:length(possible_triplets)){
            po_tr <- c(po_tr, paste(possible_triplets[[i]],collapse = ""))
          }
          possible_as <- ttoas[,po_tr]
          mut_possibilitys <- pam[ttoas[, paste(mut_triplet, collapse = "")], possible_as]
          mutation <- sample(possible_triplets, 1, prob = mut_possibilitys)
          
          if(paste(mutation[[1]], collapse = "") == paste(mut_triplet, collapse = "")){
            out <- list(sequence_complete, 0)
          } else {
            sequence_mutated <- sequence_complete
            sequence_mutated[(mut_position - pos_nucleotid + 1) : (mut_position + 3 - pos_nucleotid)] <- mutation[[1]]  
            out <- list(sequence_mutated, 1)
          }
        } else{
          out <- list(sequence_complete, 0)
        }
      }
      return(out)
    }
    
    
    
    find_true_target <- function(the_outlist, the_target_grenze){
      # Input: the_outlist: outlist mit allen RNA's und den Target informationen
      #        the_target_grenze: Grenze ab der eine mRNA ein Target wird
      
      # Output: gibt die Outliste wieder zurück, nur das in outlist[[4]] die Positionsinformationen aller Targets stehen
      
      E <- multipleE(the_outlist[[1]], the_outlist[[2]])
      trus <- which(E <= the_target_grenze)
      the_outlist[[4]]<-trus
      return(the_outlist)
    }  
    
    multipleE <- function(srna, trna) {
      # Input: srna: Eine sRNA
      #        tRNA: eine liste mit RNA's 
      # Callt IntraRNA mit gegebenen mRNA's und  einer sRNA und gibt eine liste der E-Values wieder
      # Achtung: Funktion fastaout muss im Global vorhanden sein
      # Achtung: Variable accessibillity muss im Global vorhanden sein (TRUE/FALSE)
      out <- vector("list", length(trna))
      out[] <- 0
      tempf <- tempfile()
      sreplacer <- unlist(srna)
      sreplacer <- paste(sreplacer, collapse = "")
      fasta <- fastaout(trna, tempf)
      if(accessibility){
        a <- ""
      } else {
        a <-" --qAcc=N --tAcc=N"
      }
      input <-
        paste(
          "IntaRNA --outCsvCols=id1,E -t ",
          fasta,
          " -q ",
          sreplacer,
          a,
          " --seedBP=6 --qIntLoopMax=4 --tIntLoopMax=4 --qIntLenMax=15 --tIntLenMax=15 --mode=M --outMode=C",
          collapse = "",
          sep = ""
        )
      intarna <- system(input, intern = T)
      if (length(intarna) > 1) {
        temp <- strsplit(intarna, ";")
        t1 <- lapply(temp, `[[`, 1)
        posE1 <- t1[-1]                 #Position der E-werte
        t2 <- lapply(temp, `[[`, 2)
        E1 <- t2[-1]        #E-WErte
        posE1 <- as.numeric(posE1)
        out[posE1] <- as.numeric(E1)
      }
      unlink(tempf)
      out
      
    }
    
    fastaout <- function(nlist, tempf) {
      #Input: Nucleotidliste und name des files 
      #Kreiert ein Temporaeres Fasta file aus der Nucleotidliste
      #Output: Names des Fasta files
      tempv <- c()
      sym <- ">"
      for (i in 1:length(nlist)) {
        tempc <- nlist[[i]]
        tempv2 <- paste(sym, i, collapse = "")
        tempv1 <- paste(tempc, collapse = "")
        tempv <- c(tempv, tempv2, tempv1)
      }
      a <- round(runif(1) * 10000, 0)
      writeLines(tempv, con = tempf)
      b <- tempf
      b
    }
    
    
    
    complete_value_check <- function(outlist1, mutant_list1, temp_mutation_list1, srna_value_loss1, target_value_loss1, target_grenze1, count1){
      #input: 
      #outlist and mutant_list. Outlist carrys the informations before the mutations where inserted 
      #and will also be returned if the Mutations are rejected.
      #
      #temp_mutation_list:
      
      srna <- outlist1[[1]]
      old_targets <- outlist1[[4]]
      old_srna_value <- get_value(outlist1[[2]], pos_info = old_targets)
      
      # checks if the sRNA mutated (name of sRNA = "0")
      if(names(temp_mutation_list1)[1] == "0"){
        #calculates all targets
        mutant_list1 <- find_true_target(mutant_list1, target_grenze1)
        new_target_set <- mutant_list1[[4]]
        
        lost_targets <- get_lost_targets(new_target_set, old_targets, c(1:length(outlist1[[2]])))
        #Check1: checks if there is a Value loss only looking at the old targegts
        #Check2: checks if there is a  Value loss considering all gained and lost targets
        if(target_constrain){
          check1 <- check_target_value_loss(lost_targets, old_targets, outlist1[[2]], target_value_loss1)
        
        } else {
          check1 <- 1
        }
        check2 <- check_srna_value_loss( old_srna_value, new_target_set, srna_value_loss1, outlist1[[2]])
        
        #just mrnas mutated
      } else {
        mutant_positions <- as.numeric(names(temp_mutation_list1))
        new_targets <- get_targets_tempmutlist(srna, temp_mutation_list1, target_grenze1)
        lost_targets <- get_lost_targets(new_targets, old_targets, mutant_positions)
        new_target_set <- get_new_target_set(old_targets, new_targets, lost_targets)

        
        if(target_constrain){
          check1 <- check_target_value_loss(lost_targets, old_targets, outlist1[[2]], target_value_loss1)
        } else {
          check1<- 1
        }
        check2 <- check_srna_value_loss(old_srna_value, new_target_set, srna_value_loss1, outlist1[[2]])
      }
      
      if(check1 == 1 & check2 == 1){
        mutant_list1[[4]] <- new_target_set
        mutant_list1[[3]] <- mutant_list1[[3]] + count1
        return(mutant_list1)
      } else {
        #print("no mutant")
        return(outlist1)
      }
      
      
    }
    
    get_lost_targets <- function(new_targets1, old_targets1, mutant_positions1){
      
      old_mutated_targets <- intersect(old_targets1, mutant_positions1)
      if(length(old_mutated_targets) == 0){
        return(c())
      }
      cat("Old Mutated Targets: ", old_mutated_targets)
      lost_targets <- setdiff(old_mutated_targets, new_targets1)
      return(lost_targets)
    }
    
    check_target_value_loss <- function(lost_targets1, old_targets, mrna_list, target_value_loss1){
      target_value_loss2 <- 1 - target_value_loss1
      
      if(length(lost_targets1) == 0){
        return(1)
      } 
      
      old_value <- get_value(mrna_list, pos_info =  old_targets)
      lost_value <- get_value(mrna_list, pos_info = lost_targets1)
      value_threshold <- old_value - abs(old_value * target_value_loss2)
      cat("\n lost_value: ", lost_value)
      cat("\n value_threshold: ", value_threshold)
      if((old_value - lost_value) < value_threshold){
        return(0)
      } else{
        return(1)
      }
    }
    
    get_value <-function(mrna_list1, pos_info= list(), complete_mrna_list = TRUE){
      
      value <- 0
      len <- length(pos_info)
      if(is.null(pos_info)){
        return(0)
      }
      if(length(pos_info) == 0){
        return(0)
      }
      
      for(i in 1:len){
        
        if(complete_mrna_list){
          tmp <- attr(mrna_list1[[pos_info[i]]], "value")
          value <- value + tmp
        } else{
          tmp <- attr(mrna_list1[[i]], "value")
          value <- value + tmp
        }
      }
      
      return(value)
    }
    
    get_targets_tempmutlist <- function(srna2, temp_mutation_list2, target_grenze2){
      pos <- as.numeric(names(temp_mutation_list2))
      tE <- multipleE(srna2, temp_mutation_list2)
      all_targets <- c()
      for(i in 1:length(tE)){   #Ermittlung aller Targets im Mutationsset
        tE[[i]] <- as.numeric(tE[[i]])
        if(tE[[i]] <= target_grenze2){
          all_targets <- c(all_targets, pos[i])
        }
      }
      return(all_targets)
    }
    
    get_new_target_set <- function(old_targets1, new_targets1, lost_targets1){
      temp_tar <- union(old_targets1, new_targets1)
      tar <- setdiff(temp_tar, lost_targets1)
      return(tar)
    }
    
    check_srna_value_loss <- function(old_srna_value1, new_target_set1, srna_value_loss2, mrna_list){
      new_srna_value <- get_value(mrna_list, pos_info = new_target_set1)
      if(length(new_srna_value) == 0){
        new_srna_value <- 0
      }
      #cat("\n new_srna_value: ", new_srna_value)
      value_threshold <- old_srna_value1 - abs(old_srna_value1 * (1-srna_value_loss2))
      #cat("\n srna_value_threshold: ", value_threshold)
      if(new_srna_value < value_threshold){
        return(0)
      } else{
        return(1)
      }
    }
    
    ######## Functions Ende ########
    
    ######## Initial Funktions ########
    initial_pam <- function(number = 20){
      #PAM MAtrix wird erstellt
      #Achtung: Benötig eine eine PAM1 matrix mit dem NAmen PAM1.r im Working Direktory
      
      as <-
        c(
          "ala",
          "arg",
          "asn",
          "asp",
          "cys",
          "gln",
          "glu",
          "gly",
          "his",
          "ile",
          "leu",
          "lys",
          "met",
          "phe",
          "pro",
          "ser",
          "thr",
          "trp",
          "tyr",
          "val",
          "stop"
        )
      a <- read.table("PAM1.r", header = T)
      a <- as.matrix(a)
      a <- a / 100
      
      b <- a
      for (i in 1:(number - 1)) {
        b <- b %*% a
        b <- b / 100
        
      }
      x <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      y <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100)
      b <- rbind(b, x)
      b <- cbind(b, y)
      
      pam <- b
      rownames(pam) <- as
      colnames(pam) <- as  
      return(pam)
    }
    
    make_pssm_utr <- function(sequences_complete, first_coding_nuc, pseudocount = 0){
      #sequences_complete: komplette mRNA Sequencen
      #first_coding_nuc: Der Startpunkt der Translation, Position des ersten Nucleotides nach dem Startcoden
      #
      #Erstellt eine PSSM aus allen gegebenen mRNAs. Die PSSM wird nur für den 5' UTR erstellt.
      sequences <- list()
      
      for(i in 1:length(sequences_complete)){
        sequences[[length(sequences)+1]] <- sequences_complete[[i]][0:(first_coding_nuc-1)]
      }
      mat<-matrix(,length(sequences),length(sequences[[1]]))
      for(i in 1:length(sequences)){
        mat[i,]<-tolower(sequences[[i]])  # Groß- und Kleinschreibung in den Sequenzen anpassen
      }
      
      bu<-unique(as.vector(mat))
      
      # Count Matrix erstellen
      pssm<-matrix(,length(bu),length(sequences[[1]]))
      #colnames(pssm)<-1:ncol(pssm)
      rownames(pssm)<-sort(bu)
      pssm[]<-0
      for(i in 1:ncol(mat)){
        count<-table(mat[,i])
        pssm[names(count),i]<-count
      }
      for(i in 1:ncol(pssm)){
        pssm[,i] <- ((pssm[,i] + pseudocount)/(sum(pssm[,i])+pseudocount*length(bu)))
      }
      return(pssm)
    }
    
    initial_triplet_to_as <- function(){
      #Erstellt eine Mtarix mit Hilfe derer man AS Triplets in die Aminosäuren übersetzen kann
      for (f in 1:1) {
        triplets <-
          c(
            "ttt", 
            "ttc",
            "tta",
            "ttg",
            "tct",
            "tcg",
            "tcc",
            "tca",
            "tat",
            "tac",
            "taa",
            "tag",
            "tgt",
            "tgc",
            "tga",
            "tgg",
            "ctt",
            "ctc",
            "cta",
            "ctg",
            "cct",
            "ccc",
            "cca",
            "ccg",
            "cat",
            "cac",
            "caa",
            "cag",
            "cgt",
            "cgc",
            "cga",
            "cgg",
            "att",
            "atc",
            "ata",
            "atg",
            "act",
            "acc",
            "aca",
            "acg",
            "aat",
            "aac",
            "aaa",
            "aag",
            "agt",
            "agc",
            "aga",
            "agg",
            "gtt",
            "gtc",
            "gta",
            "gtg",
            "gct",
            "gcc",
            "gca",
            "gcg",
            "gat",
            "gac",
            "gaa",
            "gag",
            "ggt",
            "ggc",
            "gga",
            "ggg"
          )
        aminosaeuren <-
          c(
            "phe",
            "phe",
            "leu",
            "leu",
            "ser",
            "ser",
            "ser",
            "ser",
            "tyr",
            "tyr",
            "stop",
            "stop",
            "cys",
            "cys",
            "stop",
            "trp",
            "leu",
            "leu",
            "leu",
            "leu",
            "pro",
            "pro",
            "pro",
            "pro",
            "his",
            "his",
            "gln",
            "gln",
            "arg",
            "arg",
            "arg",
            "arg",
            "ile",
            "ile",
            "ile",
            "met",
            "thr",
            "thr",
            "thr",
            "thr",
            "asn",
            "asn",
            "lys",
            "lys",
            "ser",
            "ser",
            "arg",
            "arg",
            "val",
            "val",
            "val",
            "val",
            "ala",
            "ala",
            "ala",
            "ala",
            "asp",
            "asp",
            "glu",
            "glu",
            "gly",
            "gly",
            "gly",
            "gly"
          )
        ttoas <- matrix(, nrow = 1, ncol = 64)
        colnames(ttoas) <- triplets
        ttoas[1, ] <- aminosaeuren
        return(ttoas)
        
      } 
    }
    
    names_to_position <- function(the_mrna, the_names){
      #Input Namen von gewünschten mRNAs
      #Gibt die Position der mRNAs innerhalb der mRNA list aus
      out<- c()
      n_rna <- names(mrna)
      if(length(n_rna) == 0){
        stop("Fehler! mRNAs besitzen keine Namen! Funktion names_to_position")
      }
      for(i in 1:length(the_names)){
        temp <- which(the_names[[i]] == n_rna)
        if(length(temp) > 0){
          out<-c(out, temp)
        }
      }
      return(out)
    }
    
    ######## Initial Functions Ende ########
    
    ######## Initialisationen ########
    temp_mutation_list <- list()
    store_count <- 0
    
    pam <- initial_pam(number = 20)
    ttoas <- initial_triplet_to_as()
    utr_pssm <- make_pssm_utr(mrna, translation_start)
    
    outlist <- vector("list", 7)
    outlist[[1]] <- srna
    outlist[[2]] <- mrna
    outlist[[3]] <- 0
    outlist <- find_true_target(outlist, target_grenze)
    outlist[[5]] <- srna
    outlist[[6]] <- outlist[[4]]
    
    
    ######## Initialsationen Ende ########
    
    ######## Wrapper ########
    
    for(iterations in 1:ceiling(durchlaeufe/store)){
      mutant_list <- outlist
      
      #Es werden X Mutationen in das sRNA/mRNA Set eingefügt
      temp_info <- stack_mutations(mutant_list, translation_start, store, rna_constrain)
      
      #mutant list ist die outlist, bloßmit eingefügten (unkontrollierten) Mutationen
      #temp_mutant_list: liste der mutierten RNAs
      mutant_list <- temp_info[[1]]
      temp_mutation_list <- temp_info[[2]]
      tmp_count <- temp_info[[3]]
      rm(temp_info)
      #Fallses keine Mutierten RNAs gibt, wird in den nächsten Iterationsschritt gesprungen.
      if(length(temp_mutation_list) == 0){
        next
      }
      #Kontrolliert ob eine sRNA mutiert ist
      
      outlist <- complete_value_check(outlist, mutant_list, temp_mutation_list, srna_value_threshold, target_value_threshold, target_grenze, tmp_count)
      
    }
    
    print("finito")
    
    ######## Wrapper Ende ########

    outlist[[7]] <- get_value(outlist[[2]], pos_info = outlist[[4]])
    outlist <- find_true_target(outlist, target_grenze)
  
    
    return(outlist)
  }
  
  ##### Funktionen Ende #####
  
  
  
  ##### Funktionen Initialisation #####
  loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  
  change_value2 <- function(srna, mrna,target_grenze, perc = 0.5, start_targets = 1){
    #This function searches for the Targets of a sRNA in a given mRNA set with a given Target threshold. 
    #Then it changes the Values of those mRNAs to percentage of the mean of all POSITIVE values.
    #Input: srna
    # List with mRNAs
    # Target_grenze = Threshold (z.B. -13)
    #perc = Percentage of the mean Value. Value assignet to the Start Targets dependent on the mean Value of all Positive Values
    # Start_targets:  Defines the number of Start Targets which will have a positive Value. All other Start Targets will have the Value 0#
    # If there are 4 start Stargets and "start_targets = 1", the first on will have a positive Value and the others a Value of 0
     
    
    #Function to find the targets of a sRNA using IntaRNA
    find_first_targets <- function(srna1, mrnas1, the_target_grenze){
      # Input: the_outlist: outlist mit allen RNA's und den Target informationen
      #        the_target_grenze: Grenze ab der eine mRNA ein Target wird
      
      # Output: gibt die Outliste wieder zurück, nur das in outlist[[4]] die Positionsinformationen aller Targets stehen
      
      E <- multipleE(srna1, mrnas1)
      trus <- which(E <= the_target_grenze)
      return(trus)
    }
    #function used by find_first_targets
    multipleE <- function(srna, trna) {
      # Input: srna: Eine sRNA
      #        tRNA: eine liste mit RNA's 
      # Callt IntraRNA mit gegebenen mRNA's und  einer sRNA und gibt eine liste der E-Values wieder
      # Achtung: Funktion fastaout muss im Global vorhanden sein
      # Achtung: Variable accessibillity muss im Global vorhanden sein (TRUE/FALSE)
      out <- vector("list", length(trna))
      out[] <- 0
      tempf <- tempfile()
      sreplacer <- unlist(srna)
      sreplacer <- paste(sreplacer, collapse = "")
      fasta <- fastaout(trna, tempf)
      if(accessibility){
        a <- ""
      } else {
        a <-" --qAcc=N --tAcc=N"
      }
      input <-
        paste(
          "IntaRNA --outCsvCols=id1,E -t ",
          fasta,
          " -q ",
          sreplacer,
          a,
          " --seedBP=6 --qIntLoopMax=4 --tIntLoopMax=4 --qIntLenMax=15 --tIntLenMax=15 --mode=M --outMode=C",
          collapse = "",
          sep = ""
        )
      intarna <- system(input, intern = T)
      if (length(intarna) > 1) {
        temp <- strsplit(intarna, ";")
        t1 <- lapply(temp, `[[`, 1)
        posE1 <- t1[-1]                 #Position der E-werte
        t2 <- lapply(temp, `[[`, 2)
        E1 <- t2[-1]        #E-WErte
        posE1 <- as.numeric(posE1)
        out[posE1] <- as.numeric(E1)
      }
      unlink(tempf)
      out
      
    }
    #function used by find_first_targets
    fastaout <- function(nlist, tempf) {
      #Input: Nucleotidliste und name des files 
      #Kreiert ein Temporaeres Fasta file aus der Nucleotidliste
      #Output: Names des Fasta files
      tempv <- c()
      sym <- ">"
      for (i in 1:length(nlist)) {
        tempc <- nlist[[i]]
        tempv2 <- paste(sym, i, collapse = "")
        tempv1 <- paste(tempc, collapse = "")
        tempv <- c(tempv, tempv2, tempv1)
      }
      a <- round(runif(1) * 10000, 0)
      writeLines(tempv, con = tempf)
      b <- tempf
      b
    }
    
    
    first_targets <- find_first_targets(srna, mrna, target_grenze)
    
    #Extract mRNA Values and get mean of all positive Values
    for(i in 1:length(mrna)){
      value_vec <- attr(mrna[[i]], "value")
    }
    pos_values <- value_vec[which(value_vec > 0)]
    mean_val <- mean(pos_values)
    #Attention bad coding here. (surely not only here)
    #Evaluates the Value the first sRNA targets should have based on the mean of the Positive Value and a given Percentage
    #If Percentage = 10 it will just take the highest positive Value 
    if(perc != 10){
      set_val <- mean_val * perc
    } else{
      set_val <- max(value_vec)
    }
    
    count<- 1
    #Replaces the Values of the sRNA Targets with the evaluated values
    for(i in first_targets){
      if(count <= start_targets){
        attr(mrna[[i]], "value") <- set_val
        count <- count + 1
      } else{
        attr(mrna[[i]], "value") <- 0
      }
    }

    return(mrna)
  }
  
  shift_value <- function(mrna, pos_border){
    #Shifts all values from a mRNA set
    #pos_border: Percentage of Positive Values in the mRNA_set
    #If pos_border = 0.3, the Values will be shifte until 30% of them are Positive (doesn`t change anything in the Value distribution)
    
    #Extract the Values from the mRNA set
    val_vec <- c()
    for(i in 1:length(mrna)){
      val_vec <- c(val_vec, attr(mrna[[i]], "value")) 
    }
    
    #Percentage of Posivtive Values
    percentage_pos <- length(which(val_vec > 0))/ length(val_vec)
    
    #Shifting all Values step by step (there is probably a smarter solution but i was too lazy to find one)
    loop <- TRUE
    while(loop){
      if(percentage_pos <= pos_border){
        loop <- FALSE
      } else {
        val_vec <- val_vec - 0.001
        percentage_pos <- length(which(val_vec > 0))/ length(val_vec)
      }
    }
    
    loop <- TRUE
    while(loop){
      if(percentage_pos >= pos_border){
        loop <- FALSE
      } else {
        val_vec <- val_vec + 0.001
        percentage_pos <- length(which(val_vec > 0))/ length(val_vec)
      }
    }
    
    #Change mRNA values to new ones
    for(i in 1:length(mrna)){
      attr(mrna[[i]], "value") <- val_vec[i]
    }
    
    return(mrna)
  }
  
  
  
  
  
  
  ##### Funktionen Initialisation Ende #####
  
  
  ##### Initialisation #####
  cores <- core_
  registerDoMC(cores)
  #Initialisierung der sRNA und der mRNA's

  mrna1 <- loadRData(mrna_)
  srna_ <- read.fasta(srna_)
  
  
  mrna1 <- shift_value(mrna1, percentage_of_pos_values)
  mrna1 <- change_value2(srna_, mrna1, target_grenze_, perc = change_start_values, start_targets = start_target_number)
  
  


  if(name_of_outputfile == "random"){
    print("You choose a Random Name for Your outputfile")
  }
  #Working direktory check
  if(save_directory != "current"){
    if(file.exists(save_directory) == FALSE){
      stop("Wrong Working directory!")
    }
    short_check<-getwd()        #kleine kontrolle ob das Working direktory auch passt
    setwd(save_directory)
    setwd(short_check)
  }
  
  ##### Initialisation Ende #####
  
  
  
  ##### Multiple Core #####
  srna <- srna_
  mrna <- mrna1
  
  daten_sammler<-list()
  for(w in 1:wiederholungen){
  #parallel simulation: 
    parallel_simulations <- foreach(i = 1:cores, .verbose = T)  %dopar% {
      rna_population <-list()
      for(n in 1:time_steps){  
        out1 <- srna_value_evolution(srna,
                                mrna,
                                translation_start_,
                                durchlaeufe_,
                                target_grenze = target_grenze_,
                                store = store_,
                                target_value_threshold = target_value_threshold_,
                                srna_value_threshold = srna_value_threshold_,
                                target_constrain = target_constrain_,
                                rna_constrain = rna_constrain_)
        #Passing the mrna and srna information along to te next run
        out1[[5]]<-srna_
        srna <- out1[[1]]
        mrna <- out1[[2]]
       # deleting the mrnas to avoid high RAM usage
        if(keep_all_rnas == FALSE){
          out1[[2]] <- "del"
        }
        #safe the generated data
        rna_population[[length(rna_population)+1]]<- out1
      }
      #return would beprobably better
      rna_poulation <- rna_population
    }

    daten_sammler[[length(daten_sammler)+1]]<-parallel_simulations
  }

  ##### Multiple Core  Ende ##### 
  
  #Function to sort the gernerated data according to the number of mutations that where inserted#
  #first row of the list: all datasets with "durchlaeufe_" numbers of Mutations
  #second row of the list: all datasets with 2 x "durchlaeufe_" numbers of Mutations
  #etc.
  daten_in_reihenfolge<-list()
  for(i in 1:length(daten_sammler)){
    for(j in 1:length(daten_sammler[[i]])){
      daten_in_reihenfolge[[length(daten_in_reihenfolge) + 1]] <- daten_sammler[[i]][[j]]
    }
  }

  daten_in_zeitlicher_reihenfolge <- list()
  for(i in 1:length(daten_in_reihenfolge[[1]])){
    
    daten_in_zeitlicher_reihenfolge[[length(daten_in_zeitlicher_reihenfolge)+1]] <- lapply(daten_in_reihenfolge, `[[`, i)
  }

  
  if(name_of_outputfile == "random"){
    name_of_outputfile <- sample(letters,5,replace = TRUE)
    name_of_outputfile <- paste(name_of_outputfile, " output", collapse ="")
    #cat("Your Output Name ist: ", name_of_outputfile)
  }
  
  if(save_directory != "current"){
    setwd(save_directory)
  }
  name_of_outputfile <- paste(name_of_outputfile,".RData",collapse="")
  save(daten_in_zeitlicher_reihenfolge, file = name_of_outputfile)
  return(daten_in_zeitlicher_reihenfolge)
}
