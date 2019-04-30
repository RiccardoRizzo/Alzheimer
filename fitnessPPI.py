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

static_hit = {}     # dictionary con l'elenco delle proteine taret dei miRNA
                    #   static_hit[miRNA] = [ [listaND], [listaD]]

ff2 = {}            # contributo di secondo livello. Per ogni miRNA si calcolano
                    # tutte le proteine collegate alle proteine target,
                    # che danno i seguenti contributi:
                    # A : prot. non dereg. collegate alle prot. target non dereg
                    # B : prot. non dereg. collegate alle prot. target dereg
                    # C : prot. dereg. collegate alle prot. target non dereg
                    # E: prot.  dereg. collegate alle prot. target  dereg

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

if __name__ == '__main__':
    lmiRNA = listaMiRNA()
    for m in lmiRNA:
        ff[m] = fitness(m)
