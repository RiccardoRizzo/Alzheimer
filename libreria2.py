import pandas as pd
import itertools

import numpy as np


#-------------------------------------------------------------------------------
def leggeListaDereg(nomeFile):
    """
    Legge il file di proteine deregolate e ricava un dict
    :param nomeFile: nome del file csv con un formato del tipo
                    numeroProgressivo, nomeProteina, peso
                    dove peso e' :
                        0 se la proteina non e' deregolata
                        1 se la proteina e' deregolata (Up o down)
    :return:    dict del tipo
                out[nomeProteina] = 0 , 1
    """

    fin = open(nomeFile, "r")
    buf= fin.readlines()
    fin.close()

    # scarta la prima riga
    del buf[0]

    out={}
    for b1 in buf:
        b=b1.strip()
        if len(b) > 0:
            try:
                nomeProte=b.strip().split(",")[0]
                valore = int( b.strip().split(",")[1] )
                out[nomeProte] = valore
            except:
                print("proteina ---", b1)
    return out


#-------------------------------------------------------------------------------
def leggeDati(nomeFilemiRNAvsProteine, nomeFileProteineDeregolate):
    """
    A partire dal file csv <nomeFilemiRNAvsProteine> con la struttura seguente:

         ,    prot1,    prot2, ... , protn
      miRNA_1,  0,      1,            1
      miRNA_2,
      ...
      miRNA_m,
      prot1,    0,       1,           0
      prot2,    1,       0 ,          0
      ...
      protm,    0,       1,           0

      dove le prime righe riportano la interazione miRNA-prot,
      e le successive riportano la interazione prot-prot
      (spazi aggiunti per chiarezza)

    e dal file <nomeFileProteineDeregolate> delle proteine,
    letto dalla funzione leggeListaDereg,

    ricava la struttura dati dd utile per il programma,
    insieme alle liste di protene e miRNA.

    :param nomeFilemiRNAvsProteine:
    :param nomeFileProteineDeregolate:

    :return:    l_p_deregolate - lista proteine deregolate
                l_p_non_deregolate - lista prot. non deregolate
                l_p_totale -  lista prot. totale
                l_miRNA_totale - lista totale dei miRNA

    """

    fin = open(nomeFileProteineDeregolate, "r")
    buf = fin.readlines()
    fin.close()

    # levo l'intestazione dal file
    buf = buf[1:]

    l_p_deregolate = []
    l_p_non_deregolate = []

    for l in buf:
        nome, d_nd = l.strip().split(",")
        if d_nd == "0":
            l_p_non_deregolate.append(nome)
        else:
            l_p_deregolate.append(nome)

    l_p_totale = l_p_deregolate + l_p_non_deregolate

    df = pd.read_csv(nomeFilemiRNAvsProteine, index_col=0)

    etichette_righe = df.index
    etichette_colonne = list(df)
    l_miRNA_totale = [ mr for mr in etichette_righe if mr not in etichette_colonne]


    return l_p_deregolate, l_p_non_deregolate, l_p_totale,  l_miRNA_totale




################################################################################

# dato un miRNA o una proteina ricava la lista delle proteine colpite
# In output du liste: proteine Deregolate (D) e non Deregolate (ND)
#-------------------------------------------------------------------------------
def hit(miRNA_Prot, df, ND, D):
    """
    Ricava la lista delle proteine collegate alla etichetta in ingresso
    input : miRNA_Prot proteina o miRNA

    output :
            ND_out : lista delle proteine Non Deregolate target dell'ingresso
            D_out : lista delle proteine deregolate target dell'ingresso

    """
    # prelevo la riga relativa al miRNA
    riga = df.loc[miRNA_Prot].values
    # lista etichette proteine
    eProt=list(df)

    lProt_out = []
    for indice in range(len(riga)):
        if riga[indice] > 0:
            lProt_out.append(eProt[indice])

    # divide in deregolate e non deregolate
    D_out = list(set(lProt_out) & set(D))
    ND_out = list(set(lProt_out) & set(ND))
    return ND_out, D_out

#-------------------------------------------------------------------------------
# calcola le proteine target dei miRNA
def static_hit(l_miRNA, df, ND, D):
    out = {}
    for miRNA in l_miRNA:
        out[miRNA] = hit(miRNA, df, ND, D)
    return out



################################################################################
################################################################################
# VRSIONE DEL FITNESS RELATIVA AL FOGLIO DEL 1 APRILE 2019 IN CUI
# OGNI PROTEINA COLLEGATA ALLE PROTEINE TARGET DEL MIRNA E' CONTATA UNA volta
# SOLA (SENZA FARE SOMME, SOLO UNIONI DI INSIEMI)
def fit_miRNA2(miRNA, df, ND, D):
    """
    Fitness del miRNA calcolato secondo il foglio di appunti
    del 1 aprile 2019. In questo caso ogni connessione e'
    contata una volta sola.
    """
    # calcolo il punteggio per le proteine colpite direttamente
    ND, D = hit(miRNA, df, ND, D)

    A=set() # insieme delle proteine collegate a ND di pk
    B=set() # insieme delle proteine colegate a D di pk
    C=set() # insieme delle proteine collegate a ND d pl
    E=set() # insieme delle proeine collegate a D di pl

    for pk in ND:
        ND_pk, D_pk = hit(pk, df, ND, D)
        A = A.union(ND_pk)
        C = C.union(D_pk)

    for pl in D:
        ND_pl, D_pl = hit(pl, df, ND, D)
        B = B.union(ND_pl)
        E = E.union(D_pl)

    return A, B, C, E

#-------------------------------------------------------------------------------
def ff2(lmiRNA, df, ND, D):
    out={}
    for miRNA in lmiRNA:
        out[miRNA] = fit_miRNA2(miRNA, df, ND, D)
    return out


################################################################################
################################################################################
##### LETTURA DEI RISULTATI ED ANALISI #########################################
################################################################################
################################################################################
#
# LETTURA DI UN file di risultati fatto in questo modo:

"""


i pesi di ottimizzazione sono (-1.0, -1.0, -1.0)

parAG:
    NGEN: 1000
    POPOLAZIONE_INIZIALE: 100
    TOURNMENT_SIZE: 20

    PROB_MUTAZIONE_SINGOLO_BIT: 0.05
    PROB_MATING: 0.8
    PROB_MUTAZIONE_INDIVIDUO: 0.2

    NUM_RIPETIZIONI : 40

# output:
#  OUTPUT_FILENAME: "fileExp.yaml"

input:
  NFILE_MIRNA_VS_PROT: "/home/riccardo/Focus/esperimenti/Anno2018/Alzheimer/Esperimenti/LavoroDiProva/prot_mirna_matrix_L0_BRCA_ER-.csv"
  NFILE_PROT_DEREG: "/home/riccardo/Focus/esperimenti/Anno2018/Alzheimer/Esperimenti/LavoroDiProva/proteins_classes_BRCA_ER-.csv"






# INIZIO RISULTATI-----------------------------------------------------------
# fine processo: tempo impiegato 0:00:16
# -----------------------------------------------------------
# fitness migliore individuo (9.0, 2466.0, 33.0)
# ELENCO miRNA-----------------------------------------------------------
# > hsa-miR-223-3p
# > hsa-miR-222-3p
# > hsa-miR-21-5p
# > hsa-miR-200b-3p
# > hsa-miR-196a-5p
# > hsa-miR-16-5p
# > hsa-miR-155-5p
# > hsa-miR-1296-5p
# > hsa-miR-20b-5p
# > hsa-miR-335-5p
# > hsa-miR-338-5p
# > hsa-miR-605-5p
# > hsa-miR-1-3p
# > hsa-miR-30a-3p
# > hsa-miR-133b
# > hsa-miR-124-3p
# > hsa-miR-30d-5p
# > hsa-miR-197-3p
# > hsa-miR-140-5p
# > hsa-miR-424-3p
# > hsa-miR-125b-5p
# > hsa-miR-183-5p
# > hsa-miR-192-5p
# > hsa-miR-92b-3p
# ELENCO PROTEINE NON COPERTE-----------------------------------------------------------
# * P02746
# * O75947
# * Q9Y2J8
# * O14579
# * P24539
# * Q16719
# * Q14956
# * P26447
# * Q9NZ01

...
seguomo altri blocchi del tipo:

# INIZIO RISULTATI-----------------------------------------------------------
...
...


"""

# i risultati sono nel file yaml di parametri, inseriti come commenti.
#
#


#-------------------------------------------------------------------------------
def cercaProva(indice, buf):
    """
    Cerca la string di inizio di una prova
    :param indice: puntatore alla lista di stringhe
    :param buf: lista di stringhe
    :return: puntatore all'inizio della lista di stringhe relativa alle prove
    """
    INIZIO_PROVA = "INIZIO RISULTATI"

    if indice >= 0:
        # cerca l'inizio della prova
        for i in range(indice, len(buf)):
            if buf[i].find(INIZIO_PROVA) > 0:
                # e' una stringa utile
                break
    else :
        i =-1
    return i





#-------------------------------------------------------------------------------
def estraeLista_miRNA(indice, buf):
    """
    Estrae la lista dei miRNA in buf a partire da indicie
    :param indice: puntatore alla lista di stringhe
    :param buf: lista di stringhe
    :return: out - lista dei miRNA
             puntatore - posizione sulla stringa successiva all'ultima letta
    """

    INIZIO_LISTA_MIRNA = "ELENCO miRNA"
    # cerca l'inizio della lista
    for i in range(indice, len(buf)):
        if buf[i].find(INIZIO_LISTA_MIRNA) > 0:
            # e' una stringa utile
            break

    out=[]
    dentro = True
    for puntatore in range(i+1, len(buf)):
        if (buf[puntatore].find(">") > 0) and dentro:
            rr= buf[puntatore].strip().split()
            out.append(rr[2])
        else :
           break
    else :
        puntatore = -1

    return out, puntatore




#-------------------------------------------------------------------------------
def estraeLista_proteine(indice, buf):
    """
    Estrae la lista delle proteine in buf a partire da indicie
    :param indice: puntatore alla lista di stringhe
    :param buf: lista di stringhe
    :return: out - lista delle proteine
             puntatore - posizione sulla stringa successiva all'ultima letta
    :return:
    """
    INIZIO_LISTA_PROTEINE = "ELENCO PROTEINE NON COPERTE"
    for i in range(indice, len(buf)):
        if buf[i].find(INIZIO_LISTA_PROTEINE) > 0:
            # e' una stringa utile
            break
    out=[]
    dentro = True
    for puntatore in range(i+1, len(buf)):
        if (buf[puntatore].find("*") > 0) and dentro:
            rr= buf[puntatore].strip().split()
            out.append(rr[2])
        else :
           break

    else:
        puntatore =-1

    return out, puntatore



#-------------------------------------------------------------------------------
def istogramma(lista):

    out={}

    for l in lista :
        if l in out.keys():
            out[l] += 1
        else :
            out[l] = 1

    return out





#-------------------------------------------------------------------------------
def filtraProteineDeregolate(dict_proteine):
    """
    Filtra il dict in uscita da leggeListaDereg(nomeFile)
     ed estrae le proteine deregolate

    :param dict_proteine:
    :return: lista di proteine deregolate
    """
    # ricavo la lista delle proteine deregolate
    l_p = []
    for kk in dict_proteine.keys():
        if dict_proteine[kk] > 0:
            l_p.append(kk)

    return l_p


#-------------------------------------------------------------------------------
def ordinaListaProteine(nomeFilemiRNAvsProteine, nomeFileProteineDeregolate):

    dd, \
    l_p_deregolate, \
    l_p_non_deregolate, \
    l_p_totale, \
    l_miRNA_utili, \
    l_miRNA_totale =    leggeDati(nomeFilemiRNAvsProteine, nomeFileProteineDeregolate)

    num_proteine_deregolate = len(l_p_deregolate)

    num_proteine_nonderegolate = len(l_p_non_deregolate)

    return l_p_totale, num_proteine_deregolate, num_proteine_nonderegolate


#-------------------------------------------------------------------------------
def creaListaOrdinataMiRNA(nomeFilemiRNAvsProteine, nomeFileProteineDeregolate):

    dd, \
    l_p_deregolate, \
    l_p_non_deregolate, \
    l_p_totale, \
    l_miRNA_utili, \
    l_miRNA_totale =    leggeDati(nomeFilemiRNAvsProteine, nomeFileProteineDeregolate)

    num_miRNA_coprenti = len(l_miRNA_utili)

    num_altri_miRNA = len(l_miRNA_totale) - len(l_miRNA_utili)

    return l_miRNA_totale, num_miRNA_coprenti, num_altri_miRNA



#-------------------------------------------------------------------------------
def miRNAvsProteine(nomeFile_miRNAvsProteine, l_miRNA, nomeFileProteineDeregolate):
    """
    Per ogni miRNA ricava la lista di poteine target divise in
    deregolate e non deregolate (target del miRNA anch'esse)

    :param nomeFile_miRNAvsProteine: file che contiene ela matrice dei dati di
                                    ingresso vedi funzione leggeDati
    :param l_miRNA: lista di miRNA risultato
    :param nomeFileProteineDeregolate: lista delle proteine deregolate
                                        (vedi funzione leggeDati)
    :return: dict del tipo
             out[<miRNA>] = [ [<lista deregolate>], [<lista non deregolate>] ]
    """

    dd, \
    l_p_tot_deregolate, \
    l_p_tot_non_deregolate, \
    l_p_totale, \
    l_miRNA_utili, \
    l_miRNA_totale =    leggeDati(nomeFile_miRNAvsProteine, nomeFileProteineDeregolate)


    # dalla lista di miRNA in ingresso ricavo quali sono le proteine su cui "impattano"
    out={}
    for miRNA in l_miRNA:
        l_etichetta = dd.loc[dd[miRNA] == 1].index.tolist()
        # ordina la lista
        l_deregolate=[]
        l_non_deregolate=[]
        for ll in l_etichetta:
            if ll in l_p_tot_deregolate:
                l_deregolate.append(ll)
            else:
                l_non_deregolate.append(ll)

        out[miRNA] = [ l_deregolate, l_non_deregolate ]

    return out


#-------------------------------------------------------------------------------
def carica_l_miRNA_risultato(nmeFilemiRNA):
    fin = open(nmeFilemiRNA, "r")
    out=fin.readlines()
    fin.close()
    # filtra tutte le strighe di lunghezza maggior edi 2
    # per evitare le righe vuote
    out=[x.strip() for x in out if len(x) > 2]
    return out
