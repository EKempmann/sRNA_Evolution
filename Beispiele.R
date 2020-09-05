#Beispiel Tree_simulation

#Es soll ein phylogenetischer Baum erstellt werden, bei dem von einem "Organismus" aus 6 Generationen generiert werden
#   wiederholungen_ = 6
#jede Generation unterscheiden 100 000 Mutationen von ihrem Vorgänger
#   durchlaeufe_ = 100000
#Jeder Organismus hat genau 2 Nachkommen
#    core_ = 2
#Als GO-Term soll iron ion binding gesetzt werden
#   go_terms_ = list(":0005506~iron ion binding,")
#(Die Angabe der GO-Terms mit":....," ist so Komplieziert, da mir Regular Expressions zu dem Zeipunkt de Programmierens noch nicht bekannt waren)
#Die Mutationen sollen nur durch die Anzahl an Targets (70% erhalt) und die PSSM/PAM20 Matrix beschränkt werden
#   term_schwelle_ = 0 UND/ODER go_term_constrain_ = FASLE, target_constrain_ = TRUE UND targets_schwelle_ = 0.7, rna_constrain_ = TRUE
#Die Grenze für ein Target soll bei -13.5 liegen
#   target_grenze_ = -13.5
#Die Räumliche Struktur der RNAs soll bei deren Interaktion mit einbezogen werden
#    accessibility = TRUE
#Es sollen 200 Mutationen angehäuft werden, bevor sie kontrolliert werden
#   store_ = 200
#DIe Output Datei soll "xyz" heißen
#   name_of_outputfile = "xyz"
#Der Output soll an einem bestimmten Ort gespeichert werden
#   safe_directory = "Pfad" (Falls es im derzeitiggen Working directory gepeichter werden soll, dann "current")
# Es sollen bestimmte mRNAs und sRNAs verwendet werden
#   mrna_ = "Datei name"
#   srna_ = "Datei name"
# Falls zufällige sRNAs oder mRNAs verwendet werden sollen, dann "random"
# Achtung: bei zufälligen mRNAs läss sich kein GO-Term Constrain einrichten
#Die Anzahl und Länge der zufälligen RNAs lässt sich mit folgenden Parametern einstellen:
#   mRNA_laenge = 40
#   mRNA_anzahl = 4357
#   sRNA_laenge = 90
#Falls die mRNA Stücke aus UTR, Starcodon und kodierendem Bereich bestehen muss der Translationsstart angegeben werden
#(Position des Ersten Nukleotides nach Startkodon)
#   translation_start_ = xx
#Falls nur UTR: translation_start_  = länge RNA, falls nur kodierender Bereich: translation_start_ = 1.
#Eine (oder mehrere)  bestimmte mRNA darf unter keinen umständen als Target verloren gehen
#   keep_mrna_ = list("b2727")
#Ansonsten
#   keep_mrna_ = "no"

go3 <- list(":0005506~iron ion binding,") 
t<- tree_sim(name_of_outputfile = "Beispiel_tree",
             save_directory = "current",
             durchlaeufe_ = 100000,            
             mrna_ = "mRNA.fa",
             srna_ = "Ryhblong.fa",
             translation_start_ = 34,
             core_ = 2,
             wiederholungen = 6,
             target_grenze_ =-13.5,
             targets_schwelle_ =0.7,
             store_ = 200,
             go_terms_ = go3,
             term_schwelle_ = 0,
             keep_mrna_ = "no",
             target_constrain_ = TRUE,
             go_term_constrain_ = FALSE,
             rna_constrain_ = TRUE,
             accessibility = TRUE,
             mRNA_laenge = 40,
             mRNA_anzahl = 4357,
             sRNA_laenge = 90)







#Beispiel Stepwise_simulation

#Es sollen 20 Durchläufe parallel gestartet werden, und in die sRNA jedes Durchlaufes sollen 10 x 25 000 Mutationen eingefügt werden.
#   time_steps = 10 
#   durchlaeufe_ = 25000
#   Die Anzahl an parallelen durchläufen lässt sich über 2 Parameter Regulieren:
#   core_ = 4
#   wiederholungen = 5
#   core gibt an wie viele Programmdurchläufe parallel auf x Kerne verteilt werden
#   wiederholungen gibt an wie oftder Prozess wiederholt werden soll. Multipliziert man beide Faktoren erhält man die Anzahl der, 
#   am Ende rauskommenden, Durchläufe.
#   Achtung: wiederholt man das ganze nur ein mal und benutzt 20 Cores, dann läuft das Programm zwar schneller, verbraucht aber viel RAM.
#Als GO-Term soll iron ion binding gesetzt werden
#   go_terms_ = list(":0005506~iron ion binding,")
#Die Mutationen sollen nur durch die Anzahl an Targets (70% erhalt) und die PSSM/PAM20 Matrix beschränkt werden
#   term_schwelle_ = 0 UND/ODER go_term_constrain_ = FASLE, target_constrain_ = TRUE UND targets_schwelle_ = 0.7, rna_constrain_ = TRUE
#Die Grenze für ein Target soll bei -13.5 liegen
#   target_grenze_ = -13.5
#Die Räumliche Struktur der RNAs soll bei deren Interaktion mit einbezogen werden
#    accessibility = TRUE
#Es sollen 200 Mutationen angehäuft werden, bevor sie kontrolliert werden
#   store_ = 200
#DIe Output Datei soll "xyz" heißen
#   name_of_outputfile = "xyz"
#Der Output soll an einem bestimmten Ort gespeichert werden
#   safe_directory = "Pfad" (Falls es im derzeitiggen Working directory gepeichter werden soll, dann "current")
# Es sollen bestimmte mRNAs und sRNAs verwendet werden
#   mrna_ = "Datei name"
#   srna_ = "Datei name"
# Falls zufällige sRNAs oder mRNAs verwendet werden sollen, dann "random"
# Achtung: bei zufälligen mRNAs lässt sich kein GO-Term Constrain einrichten
#Die Anzahl und Länge der zufälligen RNAs lässt sich mit folgenden Parametern einstellen:
#   mRNA_laenge = 40
#   mRNA_anzahl = 4357
#   sRNA_laenge = 90
#Falls die mRNA Stücke aus UTR, Starcodon und kodierendem BEreich bestehen muss der Translationsstart angegeben werden
#(Position des Ersten Nukleotides nach Startkodon)
#   translation_start_ = xx
#Falls nur UTR: translation_start_  = länge RNA, falls nur kodierender Bereich: translation_start_ = 1.
#Eine (oder mehrere)  bestimmte mRNA darf unter keinen umständen als Target verloren gehen
#   keep_mrna_ = list("b2727")
#Ansonsten
#   keep_mrna_ = "no"


go3 <- list(":0005506~iron ion binding,")
t<- stepwise_simulation2(name_of_outputfile = "95_13.5_0stepwise24",
                         save_directory = "current",    # Das Directory in dem man die Datai gespeichert haben will
                         durchlaeufe_ = 50000,         # Anzahl Mutationen pro Zyklus        
                         mrna_ = "mRNA.fa",             # mRNa filename oder "random"
                         srna_ = "Ryhblong.fa",             # sRNA filename oder "random"
                         translation_start_ = 34,
                         core_ = 4,
                         time_steps = 10,               # Anzahl an zyklen
                         wiederholungen = 6,           # wiederholungen * cores ergibt die gesamtanzahl an Wiederholungen der Simulation
                         target_grenze_ = -13.5,          # Target grenze (bei -19 hat die Eisenstess sRNA 3 anfangs Targets)
                         targets_schwelle_ = 0.95,       # Prozentsatz an Targets die behalten werden müssen
                         store_ = 200,                  # Mutationen die Gesammelt werden
                         go_terms_ = go3,               # Liste mit go Terms "no" wenn keine gesucht werden sollen
                         go_term_search = go3,          # Liste mit go_terms (hier sollte das gleiche stehen) wie in der Variable darüber
                         term_schwelle_ = 0.70,            # Prozentualer anteil an Targets mit go-Terms die behalten werden müssen ( 1 = 100%)
                         keep_mrna_ = "no",             # liste mit namen von mRNA's die behalten werden müssen
                         target_constrain_ = TRUE,
                         go_term_constrain_ = FALSE,
                         rna_constrain_ = TRUE,
                         accessibility = TRUE,   #InatRNA call mit berücksichtigung der accessibility (TRUE) oder ohne(FALSE)
                         mRNA_laenge = 40,
                         mRNA_anzahl = 4357,
                         sRNA_laenge = 25)








