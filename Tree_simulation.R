tree_sim <- function(name_of_outputfile = "random",
                            save_directory = "current",
                            durchlaeufe_ = 10000, 
                            mrna_ = "mRNA.fa",
                            srna_ = "sRNA.fa",
                            translation_start_ = 34,
                            core_ = 2,
                            wiederholungen = 1,
                            target_grenze_ = -18,
                            targets_schwelle_ = 0.6,
                            store_ = 50,
                            go_terms_ = "no",
                            term_schwelle_ = 1,
                            keep_mrna_ = "no",
                            target_constrain_ = TRUE,
                            go_term_constrain_ = TRUE,
                            rna_constrain_ = TRUE,
                            accessibility = TRUE,
                            mRNA_laenge = 40,
                            mRNA_anzahl = 4357,
                            sRNA_laenge = 25){
  ##### Beschreibung #####
  # name_of_outputfile: Name des Outputfiles (ohne .RData)
  # save_directory: gibt das Working directory an in dem der Output gespeichert werden soll
  # durchlaeufe_: Anzahl der Mutationen die pro Mutationszyklus eingefügt werden
  # mrna_:  "random" für  zufaellige mRNA's oder den namen eines files im Workingdirectory in dem die mRNA's enthalten sind
  # srna_: "random" für eine zufaellige sRNA oder den namen eines files im Workingdirectory in dem die sRNA enthalten ist
  # translations_start_: Startpunkt des Ersten nucleotides nach dem Start codon (atg)
  # core_: Anzahl der Cores die Beansprucht werden
  # wiederholungen: Angabe wie oft die Prozedur wiederholt werden soll. Die Wiederholungen entsprechen der Anzahl an "Levels" des phylogenetischen Baumes
  # target_grenze_: Grenze ab der eine mRNA ein Target ist
  # targets_schwelle_: Schwellwert an Targets die bei jeder Mutationskontrolle behalten werden müssen (0.6 = 60 %)
  # store_: Mutationen die gesammelt werden, bevor sie Kontrolliert werden
  # go_terms: Liste mit go-Terms die behalten werden müssen
  # term_schwelle: Prozentueller Anteil an Targets mit den passenden goterms, welche nach einem Mutationsdurchlauf behalten werden müssen
  # keep_mrna:  eine Liste mit den Namen der mRNA's die behalten werden MÜSSEN, dies gilt erst sobald sie targets sind.
  # target_constrain_: Gibt an ob die Mutationen unter berücksichtigung der Targets verworfen weredn sollen oder nicht (TRUE/FALSE) 
  # go_term_constrain_: Gibt an ob die Mutationen unter berücksichtigung der GO-Terms verworfen weredn sollen oder nicht (TRUE/FALSE)
  # rna_constrain: TRUE/FALSE ob auf der Mutation der mRNA ein Constrain liegen soll.
  # accessibility: TRUE/FALSE ob in der IntaRNA berechnung die accessibility der RNAs mit einbezogen werden soll
  # mrna_laenge: laenge der mRNA, falls diese Randomisiert sein sollen
  # mRNA_anzahl_ Anzahl der mRNA's falls diese Randomisiert sein sollen
  # srna_laenge: laenge der sRNA, falls diese Randomisiert sein sollen
  ##### Beschreibung Ende #####
  
  
  ##### Packages #####
  require(seqinr)
  require(doMC)
  
  ##### Packages Ende #####
  
  
  
  ##### Funktionen #####
  srna_evolution2 <- function(srna,
                              mrna,
                              translation_start,
                              durchlaeufe,
                              target_grenze = -13,
                              store = 50,
                              targets_schwelle = 0.6,
                              term_schwelle = 0.8,
                              go_terms = "no",
                              keep_mrna = "no",
                              target_constrain = TRUE,
                              go_term_constrain = TRUE,
                              rna_constrain = TRUE){
    
    #Input: srna: eine sRNA in der die einzelnen nucleotide innerhalb eines Vectors in einer Liste gespeichert sind
    #       mrna: eine Liste von mRNA's die genauso gespeichert sind wie die sRNA
    #
    #       translation_start: der Startpunkt des ersten Nucleotides nach dem Start codon (atg)
    #       druchlaeufe: Die Anzahl an Mutationen die eingebracht werden sollen
    #         Achtung: die Anzahl der tatsaechlich eingebrachten funktionen ist nicht bestimmbar.
    #         Achtung: die tatsaechliche Anzahl der Muttaionen die simuliert werden ist immer ein vielfaches  der Variable  "store"
    #                 Falls "store" 50 betraegt und die Anzahl der durchlaeufe sind 220, dann werden 250 Mutationen eingefügt
    #       target_granze: Die granze ab der eine mRNA ein Target einer sRNA ist
    #       store: Die Anzahl an Mutationen die Gesammelt werden, bevor sie selektiert werden
    #       targets_schwelle: Schwelle die angibt wie viele der alten Targets die mutierte sRNA behalten MUSS.(0.8 = 80%)
    #       term_schwelle: chwelle die angibt wie viele der Alten targets, die den Terms entsprechen, die mutierte sRNA behalten MUSS(0.8 = 80%)
    #       go_terms: Liste mit go-Terms die behalten werden müssen
    #       keep_mrna: liste mit namen der mRNA's die aufjednfall als targte behalten werden sollen
    #       target_constrain: Gibt an ob die Mutationen unter berücksichtigung der Targets verworfen weredn sollen oder nicht (TRUE/FALSE) 
    #       go_term_constrain: Gibt an ob die Mutationen unter berücksichtigung der GO-Terms verworfen weredn sollen oder nicht (TRUE/FALSE)
    #       rna_constrain: TRUE/FALSE ob auf der Mutation der mRNA ein Constrain liegen soll.
    
    
    #BESCHREIBUNG OUTLIST!!
    # Die Outlist ist ein Listen ELement mit 7 eintraegen:
    # 1.) die aktuelle sRNA
    # 2.) die aktuellen mRNA's
    # 3.) leer (artefakt)
    # 4.) Die Positionsinformationen der aktuellen targets
    # 5.) Die ursprüngliche sNRA
    # 6.) Die Ursprünglichen Targets
    # 7.) Eine Table aller Terms und deren haeufigkeit unter den Targets
    
    #Die outlist kommt in den Funktionsbeschreibungen öfters vor, da sie die meisten informationen enthaelt die Diese funktionen benötigen
    
    
    
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
      for(i in 1:store){
        current_mutation_info <- onemutation(mutant_list1[[1]], mutant_list1[[2]], translation_start1, the_rna_constrain)
        
        #Kontrolle ob die Mutation eine Veraenderung hervorgerufen hat
        if(current_mutation_info[[3]] == 0){ 
          next
        }
        
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
      return(list(mutant_list1, temp_mutation_list1))
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
    
    
    extract_one_go_term <- function(the_mrna) {
      #Input: the_mrna: eine mrna 
      #Der Name der mRNA sollte in den attributen unter "name" hinterlegt sein, da über diesen Namen die Go-Terms der mRNA ermittelt werden
      #gibt die GO terms einer mRNA aus (name irreführend! gibt alle go terms aus, nicht nur einen)
      #ACHTUNG: enviroments  hash_l2g und  hash_num2term sollten im Global verfügbar sein. (siehe Funktion: go_term_hashes)
      out <- NULL
      
      mrna_name <- attr(the_mrna, "name")
      if(length(mrna_name) == 0){
        return(list())
      }
      mrna_id <- hash_l2g[[mrna_name]]
      
      if (length(mrna_id) > 0) {
        tempi <- paste(mrna_id, collapse = "")
        termtemp <- hash_num2term[[tempi]]
        
        out <- vector("list", length(termtemp))
        if (length(termtemp) > 0) {
          for (i in 1:length(termtemp)) {
            out[[i]] <- termtemp[[i]]
          }
        }
      }
      out
    }
    
    rnas_with_current_terms <- function(the_mutant_list1, current_terms) {
      #Input: the_mutant_list1: Liste mit allen RNA's und Informationen über die Targets (siehe Outlist)
      #       current_terms: GO-Terms die immoment erhalten werden müssen
      #
      #Extrahiert alle mRNA's die immoment Targets sind und den Go-Terms entsprechen
      #Output: Liste von RNA's die den Terms entsprechen
      
      new_tru <- the_mutant_list1[[4]]
      if(length(new_tru) == 0){
        return(list())
      }
      new_targets <- the_mutant_list1[[2]][new_tru]
      rnas_with_current_terms1 <- list()
      for(i in 1:length(new_targets)){
        temp_rna_term <- extract_one_go_term(new_targets[i])
        temp_tru <- FALSE
        if(length(temp_rna_term)>0){
          for(n in 1:length(temp_rna_term)){
            if(length(which(temp_rna_term[[n]] == current_terms)) > 0){
              temp_tru <- TRUE
            }
          }
        }
        if(temp_tru){
          rnas_with_current_terms1[[length(rnas_with_current_terms1) + 1]] <- new_targets[i]
        }
      }
      return(rnas_with_current_terms1)
    }
    
    
    complete_check_srna <- function(the_outlist, the_mutant_list, mutated_srna, the_target_grenze, the_targets_schwelle, the_go_term_box, the_term_schwelle, the_keep_mrna){
      # Input: the_outlist: Liste mit allen Informationen zum allen  RNA's (siehe Outlist)
      #        the_mutated_list: Das gleiche wie the_outlist nur das hier schon alle mutationen eingefügt wurden
      #        mutated_srna: Mutierte sRNA
      #        the_target_grenze: IntaRNA Grenze ab der eine mRNA ein targe wird
      #        the_targets_schwelle: Schwelle die angibt wie viele der alten targets die mutierte sRNA behalten MUSS.(0.8 = 80%)
      #        the_go_term_box: Liste deren erster eintrag die derzeitigen go_terms sind und deren zweiter Eintrag alle target RNA's die den derzeitigen go_terms entsprechen
      #        The_term_schwelle: Schwelle die angibt wie viele der Alten targets, die den Terms entsprechen, die mutierte sRNA behalten MUSS(0.8 = 80%)
      
      # Diese Funktion ermittels zuerst die targets der mutierten sNRA unt traegt diese in the_mutated_list ein
      # Danach werden 2 checks durchgeführt: 
      # 1.) Es wird Kontrolliert ob die mutierte sRNA genug Targets behlten hat.
      #     Falls das nicht der Fall ist brcht die funktion hier direkt ab! Output ist dann -> list(0,0)
      # 2.) Es wird Kontrolliert ob die mutierte sRNA genug Targets behaelt, die den aktuellen go_terms entsprechen
      #     Falls das nicht der Fall ist brcht die funktion hier direkt ab! Output ist dann -> list(0,0)
      # Falls alle checks Positiv verlaufen ist der Output eine Liste mit 2 eintraegen: 1: 1
      #                                                                                2: the_mutated_list 
      # Anmerkung: Falls man sRNA Kontrollen hinzufügen möchte ist die beachten, dass der output dieser Funktion immer eine Liste mit 2 Elementen vorsieht,
      # wobei das erste Elemet immer eine 1 (die sRNA wird übernommen) oder eine 0 (die sRNa wird verworfen) ist und das zweite Element die alte Outlist ersezt falls die sRNA übernommen wird!
      
      
      the_mutant_list <- find_true_target(the_mutant_list, the_target_grenze)
      
      target_check <- srna_targets_check(the_outlist, the_mutant_list, the_targets_schwelle)
      
      
      if(target_check == 0){
        return(list(0,0))
      }else{
        
        if(go_terms[[1]] != "no"){
          terms_check <- srna_terms_check(the_mutant_list, the_go_term_box, the_term_schwelle)
          if(terms_check == 0){
            return(list(0,0))
          } 
        }
        
        if(the_keep_mrna[[1]] != "no"){
          spez_target_check <- srna_spezific_target_check(the_outlist, the_mutant_list, the_keep_mrna) 
        } else{
          spez_target_check <- 1
        }
        
        if(spez_target_check == 0){
          return(list(0,0))
        } else{
          return(list(1,the_mutant_list))
        }
        
        
      }
      
    }
    
    srna_targets_check <- function(the_outlist1, the_mutant_list1, the_targets_schwelle1){
      #Input: Outlidt und Mutant_list bereits beschrieben
      #       the_targets_schwelle: Schwelle für den Target erhalt 
      # Werden genug Targets behalten (z.B mehr als 80%), dann gibt die Funktion eine 1 aus, sonst eine 0.
      if(target_constrain == FALSE){
        return(1)
      }
      if(length(the_outlist1[[4]]) == 0){
        #Falls alte sRNA keine Targets hatte wird die Mutation direkt erlaubt
        return(1)
      }
      
      tru_old <- the_outlist1[[4]]
      tru_new <- the_mutant_list1[[4]]
      if(length(tru_new) == 0){
        return(0)
      }
      
      new_target_overlap<- c()
      for(i in 1:length(tru_new)){
        new_target_overlap <- c(new_target_overlap,which(tru_new[i] == tru_old))
      }
      len_new_target_overlap <- length(new_target_overlap)
      len_tru_old <-length(tru_old)
      if((len_new_target_overlap/len_tru_old) >= the_targets_schwelle1){
        return(1)
      }else{
        return(0)
      }
    }
    
    srna_spezific_target_check <- function(the_outlist1, the_mutant_list1, spezific_targets){
      # Input sind die Positionsinformationen der Spezifischen Targets!!! diese sollten als Vector vorliegen
      #
      
      old_targtets <- the_outlist1[[4]]
      if(spezific_targets[[1]] == "no"){
        return(1)
      }
      if(is.numeric(spezific_targets) == FALSE){
        stop("Spezifische Targets liegen nicht im Richtigen format vor! Funktion: srna_spezific_target_check")
      }else{
        check_spez_targ_old <- c()
        check_spez_targ_new <- c()
        for(i in 1:length(spezific_targets)){
          temp_new <- length(which(the_mutant_list1[[4]] == spezific_targets[i]))
          temp_old <- length(which(old_targtets == spezific_targets[i]))
          if (temp_new != 0){
            check_spez_targ_new <- c(check_spez_targ_new,temp_new)
          }
          if (temp_old != 0){
            check_spez_targ_old <- c(check_spez_targ_old, temp_old)
          }
        }
        
        
        if(length(check_spez_targ_new) >= length(check_spez_targ_old)){
          return(1)
        } else{
          return(0)
        }
      }
    }
    
    srna_terms_check <- function(the_mutant_list1, the_go_term_box1, the_term_schwelle1){
      #Vergleicht die Anzahl der neuen und alten Target RNAs die den Terms entsprechen
      #OUTPUT ist 1(für neue werden übernommen) oder eine 0 (für neue werden nicht übernommen) 
      if(go_term_constrain == FALSE){
        return(1)
      }
      
      old_term_rnas <- the_go_term_box1[[2]]
      if(length(old_term_rnas) == 0){
        return(1)
        
      }
      
      new_term_rnas <- rnas_with_current_terms(the_mutant_list1, the_go_term_box1[[1]])
      if(length(new_term_rnas) == 0){
        return(0)
        
      }
      names_old_term_rnas<-c()
      for(i in 1:length(old_term_rnas)){
        names_old_term_rnas <- c(names_old_term_rnas, attr(old_term_rnas[[i]], "name"))
      }
      names_new_term_rnas<-c()
      for(i in 1:length(new_term_rnas)){
        names_new_term_rnas <- c(names_new_term_rnas, attr(new_term_rnas[[i]], "name"))
      }
      
      
      
      
      same_term_rnas <- c()
      for(i in 1:length(names_new_term_rnas)){
        temp <- which(names_new_term_rnas[i] == names_old_term_rnas)
        if(length(temp) > 0){
          same_term_rnas <-c(same_term_rnas,temp)
        }
        
      }
      
      new_vs_old <- length(same_term_rnas)/length(old_term_rnas)
      if(new_vs_old >= the_term_schwelle1){
        return(1)
      } else{
        return(0)
      }
    }
    
    
    single_mrna_term_check <- function(mrna1, the_go_term_box1){
      #Kontrolliert ob eine gegebene mRNA zu den Go Terms passt
      if(length(the_go_term_box1[[1]]) == 0){
        return(0)
      } else{
        out <- 0
        mrna_term <- extract_one_go_term(mrna1)
        
        
        if(length(mrna_term)>0){
          for(l in 1:length(mrna_term)){
            if(length(which(mrna_term[[l]] == the_go_term_box1[[1]])) != 0){
              out <- 1
            }
          }
        }
        return(out)
      }
    }
    
    mrna_target_check <- function(pos_of_mutants, new_targets, old_targets){
      # Kontrolliert, ob genug Targets der sRNA erhalten bleiben
      # Output ist eine liste mit 2 Einträgen. die erste Position ist eine 1 oder 0 (Mutationen werdenübernommen/nicht übernommen).
      # Der zweite Eintrag ist eine Vectormit den Positionen der Verloren gegangenen TArgets (oder nichs, falls der Vector nicht benötigt wird).
      if(target_constrain == FALSE & go_term_constrain == FALSE & keep_mrna[[1]] == "no"){  #Constrain check
        return(list(1,integer(0)))
      }
      
      if(length(old_targets) == 0){        #Falls es bisher keine Targets gab werden alle mutationen zugelassen
        return(list(1,integer(0)))
      } 
      
      check_for_target_mutation <- match(pos_of_mutants, old_targets)
      check_for_target_mutation <- check_for_target_mutation[!is.na(check_for_target_mutation)]
      old_mutated_targets <- old_targets[check_for_target_mutation]        #Kontrolle ob unter den Mutierten Sequenzen frühere targets sind
      if(length(old_mutated_targets) == 0){                                #Falls unter den Mutationen keine Sequenz davor ein target war, werden alle Mutationen zugelassen
        return(list(1,integer(0)))
      } 
      
      new_and_old_targets <- match(new_targets, old_mutated_targets)
      new_and_old_targets<- new_and_old_targets[!is.na(new_and_old_targets)]  #Kontrolle welche der targets des Mutationssets auch davor targets waren
      if(length(new_and_old_targets) == 0){                                   # Falls ssich unter den neuen Targets keins befand, das davor ein targets war, 
        check <- (length(old_targets) - length(old_mutated_targets)) / length(old_targets)  #dann wird hier gecheckt ob die Mutationen trotz targetverlust zugelassen werden
        position_of_lost_targets <- old_mutated_targets
        if(check >= targets_schwelle | target_constrain == FALSE){
          return(list(1, position_of_lost_targets))
        } else{
          return(list(0,0))
        }
      }
      
      lost_targets <- old_mutated_targets[-new_and_old_targets]
      check <- (length(old_targets) - length(lost_targets )) / length(old_targets)  #Kontrolle wie viel Targets verloren gegangen sind
      position_of_lost_targets <- lost_targets 
      if(check >= targets_schwelle | target_constrain == FALSE){
        return(list(1, position_of_lost_targets))
      } else{
        return(list(0,0))
      }
      
    }
    
    mrna_term_check <- function(term_box, old_tar1, position_of_lost_targets, all_mrna){
      #Kontrolliert ob genug der Targets, die dem GO-Term entsprechen, erhalten geblieben sind
      # Output ist eine 1 (genug sind erhalten), oder eine 0 (nicht genug erhalten)
      
      if(length(position_of_lost_targets) == 0 | go_term_constrain == FALSE){
        return(1)
      }
      check_term <- c()
      if(go_terms[[1]] != "no"){
        for(i in 1:length(old_tar1)){
          temp <- single_mrna_term_check(all_mrna[[old_tar1[i]]], term_box)
          check_term <- c(check_term, temp)
        }
      } else{
        return(1)
      }
      
      position_of_lost_targets <- match(position_of_lost_targets, old_tar1)
      position_of_lost_targets <- position_of_lost_targets[!is.na(position_of_lost_targets)]
      
      num_of_lost_go_term <- sum(check_term[position_of_lost_targets])
      num_of_rnas_with_terms <- length(term_box[[2]])
      vgl <- (num_of_rnas_with_terms - num_of_lost_go_term) /num_of_rnas_with_terms 
      if(num_of_rnas_with_terms == 0){
        return(1)
      } else if(vgl >= term_schwelle){
        return(1)
      } else {
        return(0)
      }
      
    }
    
    check_spez_traget_mrna <- function(spezific_targets, position_of_lost_targets){
      #checkt für jedes spezifische Target ob es in dem Vector der Verloren gegangenen Targets ist
      
      for(i in spezific_targets){
        if(i %in% position_of_lost_targets){
          return(0)
        }
      }
      return(1)
      
    }
    
    complete_check_mrna <- function(the_outlist1, the_mutant_list1, the_temp_mutation_list1, the_target_grenze1, the_go_term_box, keep_mrna){
      # Falls nur mRNAs mutiert sind, wird hier kontrolliert ob die Mutationen den vorgegebenen Constrains entsprechen, 
      #falls ja werden sie übernommen und die Outlist wird aktualisiert, falls nicht werden sie Verworfen.
      #####  Vorbereitung  #####
      srna1 <- the_outlist1[[1]][1]
      mrnas <- the_outlist1[[2]]
      
      #Positionen der Mutationen werden extrahiert und Bindungsenergien zur sRNA berechnet
      pos <- as.numeric(names(the_temp_mutation_list1))
      tE <- multipleE(srna1,the_temp_mutation_list1)
      old_tar <- the_outlist1[[4]]
      
      new_targets <- c()
      for(i in 1:length(tE)){   #Ermittlung aller Targets im Mutationsset
        tE[[i]] <- as.numeric(tE[[i]])
        if(tE[[i]] <= the_target_grenze1){
          new_targets <- c(new_targets, pos[i])
        }
      }
      #####  Vorbereitung Ende  #####
      
      
      
      #####  Wrapper  #####
      
      tmp <- mrna_target_check(pos, new_targets, old_tar)
      
      check_tar <- tmp[[1]]
      pos_lost_targets <- tmp[[2]]
      rm(tmp)
      if(check_tar == 0){
        return(the_outlist1)
      }
      check_term <- mrna_term_check(the_go_term_box, old_tar, pos_lost_targets, mrnas)
      if(check_term == 0){
        return(the_outlist1)
      }
      
      if(keep_mrna[[1]] != "no"){
        check_spz_target <- check_spez_traget_mrna(keep_mrna, pos_lost_targets)
        if(check_spz_target == 0){
          return(the_outlist1)
        }
      }
      
      #Generierung des Outputs
      cut_vec1 <- match(new_targets,old_tar)
      cut_vec1  <- cut_vec1[!is.na(cut_vec1)]
      if(length(cut_vec1) != 0){
        old_tar <- old_tar[-cut_vec1]
      }
      
      cut_vec2 <- match(pos_lost_targets, old_tar)
      cut_vec2  <- cut_vec2[!is.na(cut_vec2)]
      if(length(cut_vec2) != 0){
        old_tar <- old_tar[-cut_vec2]
      }
      tar <- c(old_tar, new_targets)
      
      boxcheck <<- TRUE
      the_outlist1 <- the_mutant_list1
      the_outlist1[[4]] <- tar
      #####  Wrapper Ende  #####
      return(the_outlist1)
    }
    
    
    
    
    
    ######## Functions Ende ########
    
    ######## Initial Funktions ########
    get_current_terms <- function(the_outlist,additional_terms = list(), dontdo = TRUE){
      
      
      if(dontdo){
        return(additional_terms)
      }
      out_terms <- list()
      trus <- the_outlist[[4]]
      if(length(trus) > 0){  
        for(i in 1:length(trus)){
          temp_term <- extract_one_go_term(the_outlist[[2]][[trus[i]]])
          if(length(temp_term) > 0){
            out_terms[[length(out_terms)+1]] <- temp_term 
          }
        }
      }
      out_terms <- c(additional_terms,out_terms)
      return(out_terms)
    }
    
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
    
    go_term_hashes <- function(genid_2_annotation = "genid_2_anno.txt", locustag2geneid = "locustag2geneid.txt"){
      #Hashes Vorbereiten hash_num2term:Jeder Nummer sind Terme zugeordnet  hash_term2num: Jedem Term sind Nummern zugeortnet hash_l2g: Jedem Name sind Nummern Zugeordnet 
      #Es müssen die Dataien genid_2_anno.txt und locustag2geneid.txt im Working directory liegen
      
      a<-read.csv(genid_2_annotation,sep="\t",header=T)
      b<-a[[1]]
      c1<-as.vector(a[[5]])
      c2<-as.vector(a[[6]])
      c3<-as.vector(a[[7]])
      
      hash_num2term1<-new.env(hash=T)
      hash_term2num1<-new.env(hash=T)
      
      for(i in 1:length(b)){
        tempb<-paste(b[[i]],collapse="") 
        tempc1<-strsplit(as.vector(c1[[i]]),"GO")
        tempc2<-strsplit(as.vector(c2[[i]]),"GO")
        tempc3<-strsplit(as.vector(c3[[i]]),"GO")
        hash_num2term1[[tempb]]<-c(tempc1[[1]][-1],tempc2[[1]][-1],tempc3[[1]][-1])
      }
      for(i in 1:length(b)){
        tempb<-paste(b[[i]],collapse="") 
        tempc1<-strsplit(as.vector(c1[[i]]),"GO")
        tempc2<-strsplit(as.vector(c2[[i]]),"GO")
        tempc3<-strsplit(as.vector(c3[[i]]),"GO")
        tempc1<-tempc1[[1]][-1]
        tempc2<-tempc2[[1]][-1]
        tempc3<-tempc3[[1]][-1]
        if(length(tempc1!=0)){
          for(j in 1:length(tempc1)){
            hash_term2num1[[tempc1[j]]]<-c(hash_term2num1[[tempc1[j]]],tempb)
          }
        }
        if(length(tempc2!=0)){
          for(j in 1:length(tempc2)){
            hash_term2num1[[tempc2[j]]]<-c(hash_term2num1[[tempc2[j]]],tempb)
          }
        }
        if(length(tempc3!=0)){
          for(j in 1:length(tempc3)){
            hash_term2num1[[tempc3[j]]]<-c(hash_term2num1[[tempc3[j]]],tempb)
          }
        }
      }
      
      
      locus2genid<-read.csv(locustag2geneid,sep="\t",header=T)
      l2ga<-locus2genid[[1]]
      l2gb<-locus2genid[[2]]
      hash_l2g1<-new.env(hash=T)
      for(i in 1:length(l2gb)){
        tempb<-paste(l2ga[[i]],collapse="") 
        tempc<-l2gb[[i]]
        hash_l2g1[[tempb]]<-tempc
      }
      
      return(list(hash_num2term1, hash_term2num1, hash_l2g1))
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
    
    pam <- initial_pam()
    ttoas <- initial_triplet_to_as()
    utr_pssm <- make_pssm_utr(mrna, translation_start)
    temp_hash<-go_term_hashes()
    hash_num2term <- temp_hash[[1]]
    hash_term2num <- temp_hash[[2]]
    hash_l2g <- temp_hash[[3]]
    rm(temp_hash)
    
    
    outlist <- vector("list", 7)
    outlist[[1]] <- srna
    outlist[[2]] <- mrna
    outlist <- find_true_target(outlist, target_grenze)
    outlist[[5]] <- srna
    outlist[[6]] <- outlist[[4]]
    
    current_terms <- get_current_terms(outlist, additional_terms = go_terms)
    
    go_term_box<-list()
    go_term_box[[1]] <- current_terms
    go_term_box[[2]] <- rnas_with_current_terms(outlist, current_terms)
    boxcheck <- FALSE
    if(keep_mrna[[1]] != "no"){
      keep_mrna <- names_to_position(mrna, keep_mrna)
    }
    
    
    
    
    
    
    
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
      rm(temp_info)
      #Fallses keine Mutierten RNAs gibt, wird in den nächsten Iterationsschritt gesprungen.
      if(length(temp_mutation_list) == 0){
        next
      }
      
      #Kontrolliert ob eine sRNA mutiert ist
      
      if(names(temp_mutation_list)[1] == "0"){
        #sRNA ist Mutiert
        mutierte_srna <- temp_mutation_list[[1]]
        check_srna <- complete_check_srna(outlist, mutant_list, mutierte_srna, target_grenze, targets_schwelle, go_term_box, term_schwelle, keep_mrna)
        if(check_srna[[1]] == 1){
          outlist <- check_srna[[2]]
          go_term_box[[2]] <- rnas_with_current_terms(outlist, current_terms)
        }
      } else{
        #keine sRNA ist Mutiert
        outlist <- complete_check_mrna(outlist, mutant_list, temp_mutation_list, target_grenze, go_term_box, keep_mrna)
        if(boxcheck){
          go_term_box[[2]] <- rnas_with_current_terms(outlist, current_terms)
          boxcheck <- FALSE
        }
      }
    }
    
    
    
    ######## Wrapper Ende ########
    
    
    terms_end <- get_current_terms(outlist, dontdo = FALSE)
    outlist[[7]]<- table(unlist(terms_end))
    rm(terms_end)
    
    outlist <- find_true_target(outlist, target_grenze)
    
    
    return(outlist)
  }
  
  
  
  ##### Funktionen Ende #####
  
  
  
  ##### Funktionen Initialisation #####
  permu <- function(lae, an) {
    out <- vector("list", an)
    nucleotid <- c("c","a","t","g")
    vec <- function(lae) {
      pl <- c()
      for (i in 1:lae) {
        pl <-c(pl,sample(nucleotid,1))  
      }
      return(pl)
    }
    for (i in 1:an) {
      out[[i]] <- vec(lae)
    }
    return(out)
  }
  ##### Funktionen Initialisation Ende #####
  
  
  
  ##### Initialisation #####
  cores <- core_
  registerDoMC(cores)
  #Initialisierung der sRNA aun der mRNA's
  #Falls man Random sRNA/mRNAs will, werden diese hier initialisiert
  if (mrna_ == "random") {
    mrna_ <- permu(mRNA_laenge, mRNA_anzahl)
  } else{
    mrna_ <- read.fasta(mrna_)
  }
  
  if (srna_ == "random") {
    srna_ <- permu(sRNA_laenge, 1)
  } else{
    srna_ <- read.fasta(srna_)
  }
  
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
  rna_seq <-list(list(srna_,mrna_))
  daten_sammler<- vector("list", wiederholungen)
  
  for(w in 1:wiederholungen){
    tmp <- list()
    for(m in 1:length(rna_seq)){
      parallel_simulations <- foreach(i = 1:cores, .verbose = T)  %dopar% {
        out2 <- srna_evolution2(rna_seq[[m]][[1]],
                               rna_seq[[m]][[2]],
                               translation_start_,
                               durchlaeufe_,
                               target_grenze = target_grenze_,
                               store = store_,
                               targets_schwelle = targets_schwelle_,
                               term_schwelle = term_schwelle_,
                               go_terms = go_terms_,
                               target_constrain = target_constrain_,
                               go_term_constrain = go_term_constrain_,
                               rna_constrain = rna_constrain_,
                               keep_mrna = keep_mrna_)
      }
      
      for(p in 1:length(parallel_simulations)){
        tmp[[length(tmp) + 1]] <- list(parallel_simulations[[p]][[1]], parallel_simulations[[p]][[2]])
        daten_sammler[[w]][[length(daten_sammler[[w]])+1]] <- parallel_simulations[[p]]
      }
      
    }
    rna_seq <- tmp
  }
  ##### Multiple Core  Ende ##### 
  
  
  
  
  #####
  ####
  
  
  ##
  
  if(name_of_outputfile == "random"){
    name_of_outputfile <- sample(letters,5,replace = TRUE)
    name_of_outputfile <- paste(name_of_outputfile, " output", collapse ="")
    cat("Your Output Name ist: ", name_of_outputfile)
  }
  
  if(save_directory != "current"){
    setwd(save_directory)
  }
  name_of_outputfile <- paste(name_of_outputfile,".RData",collapse="")
  save(daten_sammler, file = name_of_outputfile)
  return(daten_sammler)
}