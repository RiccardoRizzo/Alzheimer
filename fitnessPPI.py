import pandas as pd
import numpy as np

# import matplotlib.pyplot as plt; plt.rcdefaults()

#import matplotlib.pyplot as plt



########################################################################
########################################################################
### riferimenti vuoti da riempire a cura del programma
D = []              # lista proteine deregolate
ND = []             # lista proteine non deregolate
df=pd.DataFrame()   # dataframe della struttura miRNA -PPI

static_hit = {}     # dictionary con l'elenco delle proteine target dei miRNA
                    #   static_hit[miRNA] = [ [listaND], [listaD] ]

ff2 = {}            # contributo di secondo livello. Per ogni miRNA si calcolano
                    # tutte le proteine collegate alle proteine target,
                    # che danno i seguenti contributi:
                    # A : prot. non dereg. collegate alle prot. target non dereg
                    # B : prot. non dereg. collegate alle prot. target dereg
                    # C : prot. dereg. collegate alle prot. target non dereg
                    # E : prot. dereg. collegate alle prot. target  dereg

listamiRNA = []     # lista dei miRNA nel problema
########################################################################
########################################################################

#-------------------------------------------------------------------------------
def listaProtTarget(individual):
    # si uniscono gli insiemi di tutte le proteine deregolate e non Deregolate
    # e si contano gli elementi
    ND_T = set()
    D_T = set()
    for ii in range(len(individual)):
        if individual[ii] > 0:
            ND, D = static_hit[listamiRNA[ii]]
            ND_T = ND_T.union(ND)
            D_T = D_T.union(D)

    return  D_T, ND_T



#-------------------------------------------------------------------------------
def fitness_L0(individual):
    l_miRNA = listamiRNA

    # fitness di primo livello
    #----------------------------
    # si uniscono gli insiemi di tutte le proteine deregolate e non Deregolate
    # e si contano gli elementi
    D_T, ND_T = listaProtTarget(individual)

    out_D_T = len(D_T)
    out_ND_T = len(ND_T)

    return  out_D_T, out_ND_T


def fitness_L0_num(individual):
    """
    Aggiunge all fitness di livello zero anche il numero di miRNA
    di cui e' composto l'individuo (parametro da minimizzare)
    """
    D_T, ND_T =  fitness_L0(individual)
    temp_array = np.array(individual)
    lun = temp_array.sum()
    return D_T, ND_T, lun





def fitness_L0_L1(individual, R=0.5):
    """
    Calcola il fitness della lista di miRNA passata
    input: lmiRNA: lista dei miRNA di cui clacolare il Fitness
    output: out : punteggio complessivo della lista di miRNA
    """
    l_miRNA = listamiRNA

    # fitness di primo livello
    #----------------------------
    # si uniscono gli insiemi di tutte le proteine deregolate e non Deregolate
    # e si contano gli elementi
    out_D_T, out_ND_T  = fitness_L0(individual)

    # fitness di secondo livello
    #----------------------------
    # si considerano tutte le proteine collegate all'insieme precedente
    # e si esegue l'unione. quindi si contano gli elementi dell'insieme.
    A_T=set() # insieme delle proteine collegate a ND di pk
    B_T=set() # insieme delle proteine colegate a D di pk
    C_T=set() # insieme delle proteine collegate a ND d pl
    E_T=set() # insieme delle proeine collegate a D di pl

    for ii in range(len(individual)):
        if individual[ii] > 0:
            A, B, C, E =  ff2[l_miRNA[ii]]
            A_T = A_T.union(A)
            B_T = B_T.union(B)
            C_T = C_T.union(C)
            E_T = E_T.union(E)

    out_D_T += R * ( len(C_T) + len(E_T) )
    out_ND_T += R * ( len(A_T) + len(B_T) )

    return  out_D_T, out_ND_T


def fitness_L0_L1_num(individual):
    """
    Aggiunge alla fitness di livello uno anche il numero di miRNA
    di cui e' composto l'individuo (parametro da minimizzare)
    """
    D_T, ND_T =  fitness_L0_L1(individual)
    temp_array = np.array(individual)
    lun = temp_array.sum()
    return D_T, ND_T, lun








#-------------------------------------------------------------------------------
def fitness2(individual, R=0.5):
    """
    Calcola il fitness della lista di miRNA passata
    input: lmiRNA: lista dei miRNA di cui clacolare il Fitness
    output: out : punteggio complessivo della lista di miRNA
    """
    l_miRNA = listamiRNA
    out = 0

    # fitness di primo livello
    #----------------------------
    # si uniscono gli insiemi di tutte le proteine deregolate e non Deregolate
    # e si contano gli elementi
    ND_T = set()
    D_T = set()
    for ii in range(len(individual)):
        if individual[ii] > 0:
            ND , D = static_hit[l_miRNA[ii]]
            ND_T = ND_T.union(ND)
            D_T = D_T.union(D)
    out = len(D_T) - len(ND_T)

    # fitness di secondo livello
    #----------------------------
    # si considerano tutte le proteine collegate all'insieme precedente
    # e si esegue l'unione. quindi si contano gli elementi dell'insieme.
    A_T=set() # insieme delle proteine collegate a ND di pk
    B_T=set() # insieme delle proteine colegate a D di pk
    C_T=set() # insieme delle proteine collegate a ND d pl
    E_T=set() # insieme delle proeine collegate a D di pl

    for ii in range(len(individual)):
        if individual[ii] > 0:
            A, B, C, E =  ff2[l_miRNA[ii]]
            A_T = A_T.union(A)
            B_T = B_T.union(B)
            C_T = C_T.union(C)
            E_T = E_T.union(E)
    out += R*( len(C_T) - len(A_T) + len(E_T) - len(B_T) )

    return out,

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
### ==========================================================================
###  NUOVE FUNZIONI OBIETTIVO ================================================
### ==========================================================================
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

# struttura di l_miRNA 
    
    # la prima lista indica la copertura sulle proteine deregolate
    # il secondo numero indica il peso.
    #
    # l_miRNA = [
    #     [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 0],
    #     [[0, 0, 0, 0, 0, 0, 0, 0, 0, 1], 1],
    #     [[0, 0, 0, 0, 0, 0, 0, 0, 1, 0], 2],
    #     [[0, 0, 0, 0, 0, 0, 0, 0, 1, 1], 3],
    #     [[0, 0, 0, 0, 0, 0, 0, 1, 0, 0], 4],
    #     [[0, 0, 0, 0, 0, 0, 0, 1, 0, 1], 5],
    #     [[0, 0, 0, 0, 0, 0, 0, 1, 1, 0], 6],
    #     [[0, 0, 0, 0, 0, 0, 0, 1, 1, 1], 7],
    #     [[0, 0, 0, 0, 0, 0, 1, 0, 0, 0], 8],
    #     [[0, 0, 0, 0, 0, 0, 1, 0, 0, 1], 9],
    #     [[0, 0, 0, 0, 0, 0, 1, 0, 1, 0], 10]
    # ]




def num_proteine_coperte(individual):
    """
    Data la popolazione descritta in individual, 
    calcola il numero di proteine coperte

        :param individual: vettore di zero e 1 di lunghezza 
                            pari al numero totale di miRNA
        :param l_miRNA: lista di strutture con la copertura
                    sulle proteine ed il peso
    """
    D_T, ND_T = fitness_L0(individual)
    valore = len(D_T)
    return valore

#--------------------------------------------------------------------------
def num_proteine_coperte_piu_di_K_volte(individual, l_miRNA, K):
    """
    Data la popolazione descritta in individual, 
    calcola il numero di proteine coperte piu' di K volte

        :param individual: vettore di zero e 1 di lunghezza 
                            pari al numero totale di miRNA
        :param l_miRNA: lista di strutture con la copertura
                    sulle proteine ed il peso
    """
    # crea un vettore per memorizzare la copertua delle proteine
    ss=np.zeros(len(l_miRNA[0][0],))

    for ii in range(len(individual)):
        if individual[ii] > 0:
            ss = np.add(ss, l_miRNA[ii][0])

    valore = (ss > K).sum()

    return valore


#--------------------------------------------------------------------------
def peso_totale_soluzione(individual, l_miRNA):
    """
    Calcola il peso totale della soluzione
    """
    peso_totale = 0
    for ii in range(len(individual)):
        if individual[ii] > 0:
            peso_totale += l_miRNA[ii][1]
    return peso_totale



#--------------------------------------------------------------------------
def num_miRNA_coprenti_piu_di_K_prot(individual, l_miRNA, K):
    """
        Data la popolazione descritta in individual, 
        calcola il numero di miRNA 
        che coprono piu' di K proteine

        :param individual: vettore di zero e 1 
                            di lunghezza pari al numero totale di miRNA
        :param l_miRNA: lista di strutture con la copertura
                    sulle proteine ed il peso

        :param K : desiderato numero di 
                    proteine coperte dal singolo miRNA
    """
    totale = 0
    for ii in range(len(individual)):
        if individual[ii] > 0:
            kk = np.sum(l_miRNA[ii][0])
            if kk > K:
                totale += 1
    
    return totale


#--------------------------------------------------------------------------
def num_miRNA_di_peso_minore_a_K(individual, l_miRNA, K):
    """
        Data la popolazione descritta in individual, 
        calcola il numero di miRNA di peso inferiore a K

        :param individual: vettore di zero e 1 
                            di lunghezza pari al numero totale di miRNA
        :param l_miRNA: lista di strutture con la copertura
                    sulle proteine ed il peso

        :param K : peso massimo dei miRNA
    """
    totale = 0
    for ii in range(len(individual)):
        if individual[ii] > 0 and l_miRNA[ii][1] < K:
            totale += 1 
    return totale


#--------------------------------------------------------------------------
def obiettivo2(individual, l_miRNA):
    """
    Calcola il punteggio della soluzione

    :param individual: vettore di zero e 1 di lunghezza eguale al numero totale di miRNA
    :param l_miRNA: lista di strutture con la copertura
                    sulle proteine ed il peso

    :return:    valore - indica il numero di proteine coperte dalla soluzione
                peso - indica il numero di proteine non deregolate imapttate
                        dalla soluzione
                svantaggio_ss - indica il numero di proteine coperte piu' volte
                                dalla soluzione
    """
    valore =  num_proteine_coperte(individual, l_miRNA)
    peso =  peso_totale_soluzione(individual, l_miRNA)
    svantaggio_ss = num_proteine_coperte_piu_di_K_volte(individual, l_miRNA, K=1)

    return valore, peso, svantaggio_ss











##############################################################################
## MAIN  #####################################################################
##############################################################################

if __name__ == '__main__':
    #lmiRNA = listamiRNA
    # for m in lmiRNA:
    #     ff[m] = fitness(m)
    
    #####################################################################
    #####################################################################
    ####### COLLAUDO DELLA FITNESS MODULARE
    ####################################################################

    
    l_miRNA = [
        [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 0],
        [[0, 0, 0, 0, 0, 0, 0, 0, 0, 1], 1],
        [[0, 0, 0, 0, 0, 0, 0, 0, 1, 0], 2],
        [[0, 0, 0, 0, 0, 0, 0, 0, 1, 1], 3],
        [[0, 0, 0, 0, 0, 0, 0, 1, 0, 0], 4],
        [[0, 0, 0, 0, 0, 0, 0, 1, 0, 1], 5],
        [[0, 0, 0, 0, 0, 0, 0, 1, 1, 0], 6],
        [[0, 0, 0, 0, 0, 0, 0, 1, 1, 1], 7],
        [[0, 0, 0, 0, 0, 0, 1, 0, 0, 0], 8],
        [[0, 0, 0, 0, 0, 0, 1, 0, 0, 1], 9],
        [[0, 0, 0, 0, 0, 0, 1, 0, 1, 0], 10],
    ]

    individual = [1,1,1,0,0,0,0,0,0,0,1]


    print("num_proteine_coperte=", num_proteine_coperte(individual, l_miRNA))
    print("num_proteine_coperte_piu_di_K_volte=", num_proteine_coperte_piu_di_K_volte(individual, l_miRNA, 1))
    print("peso_totale_soluzione=", peso_totale_soluzione(individual, l_miRNA))
    print("num_miRNA_coprenti_piu_di_K_prot=", num_miRNA_coprenti_piu_di_K_prot(individual, l_miRNA, 1))
    print("num_miRNA_di_peso_minore_a_K=", num_miRNA_di_peso_minore_a_K(individual, l_miRNA, 3))
    print("obiettivo2=", obiettivo2(individual, l_miRNA))
