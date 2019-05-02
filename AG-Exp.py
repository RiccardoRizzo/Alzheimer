import pandas as pd
import numpy as np



import random

from deap import algorithms
from deap import base
from deap import creator
from deap import tools

import sys


# import funzione_obiettivo as fo
import fitnessPPI as fo
import libreria2 as lb

import datetime

#import libPar
#import datiEmail as dm
import yaml

import os.path

"""
In questo programma la ottimizzazione e' fatta considerando
il problema come weighted cover set.

Le chiamate a funzione e alcune variabili sono state cambiate in accordo.
"""

# ------------------------------------------------------------------------------
def cxTwoPointCopy(ind1, ind2):
    """Execute a two points crossover with copy on the input individuals. The
    copy is required because the slicing in numpy returns a view of the data,
    which leads to a self overwritting in the swap operation. It prevents
    ::

        >>> import numpy
        >>> a = numpy.array((1,2,3,4))
        >>> b = numpy.array((5.6.7.8))
        >>> a[1:3], b[1:3] = b[1:3], a[1:3]
        >>> print(a)
        [1 6 7 4]
        >>> print(b)
        [5 6 7 8]
    """
    size = len(ind1)

    cxpoint1 = random.randint(1, size)
    cxpoint2 = random.randint(1, size - 1)
    if cxpoint2 >= cxpoint1:
        cxpoint2 += 1
    else:  # Swap the two cx points
        cxpoint1, cxpoint2 = cxpoint2, cxpoint1

    ind1[cxpoint1:cxpoint2], ind2[cxpoint1:cxpoint2] \
        = ind2[cxpoint1:cxpoint2].copy(), ind1[cxpoint1:cxpoint2].copy()

    return ind1, ind2




#-------------------------------------------------------------------------------
def protDeregNonCoperte(individual, l_p_deregolate):
    D, ND = fo.listaProtTarget(individual)

    # sottrae dall'insieme l_p_deregolate l'insieme D
    return set(l_p_deregolate).difference(D)


# ------------------------------------------------------------------------------
def main(filePar):

    # carica i valori da passare ai programmi
    with open(filePar, 'r') as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)

        NGEN = cfg["parAG"]["NGEN"]
        POPOLAZIONE_INIZIALE = cfg["parAG"]["POPOLAZIONE_INIZIALE"]
        TOURNMENT_SIZE = cfg["parAG"]["TOURNMENT_SIZE"]
        PROB_MUTAZIONE_SINGOLO_BIT = cfg ["parAG"]["PROB_MUTAZIONE_SINGOLO_BIT"]
        PROB_MATING = cfg ["parAG"]["PROB_MATING"]
        PROB_MUTAZIONE_INDIVIDUO = cfg["parAG"]["PROB_MUTAZIONE_INDIVIDUO"]

        NUM_RIPETIZIONI = cfg["parAG"]["NUM_RIPETIZIONI"]

        fitness = cfg["parAG"]["FUNZIONE_FITNESS"]
        fitness_eval = cfg["parAG"]["FITNESS_EVAL"]

        nomeDataset = cfg["input"]["NFILE_MIRNA_VS_PROT"]

        nomeF_protDereg = cfg["input"]["NFILE_PROT_DEREG"]


    # random.seed(64)

    ############# PREPARAZIONE AMBIENTE fo #####################################
    fo.df = pd.read_csv(nomeDataset, index_col=0)
    # PREPARAZIONE DATASET
    l_p_deregolate, \
    l_p_non_deregolate, \
    l_p_totale, \
    l_miRNA_totale = lb.leggeDati(nomeDataset, nomeF_protDereg)

    fo.ND = l_p_non_deregolate
    fo.D = l_p_deregolate
    fo.listamiRNA = l_miRNA_totale

    fo.ff2 = lb.ff2(l_miRNA_totale, fo.df, fo.ND, fo.D)
    fo.static_hit = lb.static_hit(l_miRNA_totale, fo.df, fo.ND, fo.D)
    ############################################################################

    print("fine inizializzazione ambiente")
    print("NOME DEL DATASET", nomeDataset)
    print("numero proteine deregolate ", str(len(l_p_deregolate)) )

    # INIZIO PARTE ALGORITMO GENETICO
    #DIM_INDIVIDUO = len(l_p)
    DIM_INDIVIDUO = len(l_miRNA_totale)

    #OUTPUT_FILENAME = cfg["output"]["OUTPUT_FILENAME"]
    # al momento l'output e' scritto in coda al file parametri
    OUTPUT_FILENAME = filePar

    #
    # La minimizzazione e' sulle tre uscite della funzione obiettivo:
    # valore, - numero di zeri nella lista somma, quindi numero delle
    #           proteine non coperte
    # peso, - numero degli "impatti" sulle proteine non deregolate
    # svantaggio_ss - numero degli elementi maggiori di 1 nella lista somma,
    #                   quindi numero delle proteine che risultano coperte
    #                   piu' volte nel risultato
    # mettendo la tupla a (-1.0, -0.5, -0.01) uso solo due termini
    # creator.create("FitnessMax", base.Fitness, weights=(-1.0, -0.5, -0.0001))
    # creator.create("FitnessMax", base.Fitness, weights=(-1.0, -1.0, -1.0))

    # fitness piu' semplice per obiettivo3
    # creator.create("FitnessMax", base.Fitness, weights=(-1.0 , ))
    # fitness piu' complessa per obiettivo 2
    #creator.create("FitnessMax", base.Fitness, weights=(-1.0,))

    # se si deve massimizzare l'obiettivo
    # creator.create("FitnessMax", base.Fitness, weights=(1.0,))

    # il fitness adesso ha due elementi: massimizzare l'impatto sulle deregolate
    # minimizzare l'impatto sulle non deregolate
    #creator.create("FitnessMax", base.Fitness, weights=(1.0,-1.0))
    
    creator.create("FitnessMax", base.Fitness, weights=fitness_eval)
    creator.create("Individual", np.ndarray, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()
    toolbox.register("attr_bool", random.randint, 0, 1)
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, n=DIM_INDIVIDUO)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    ### DEFINIZIONE DELLA FUNZONE DI FITNESS -----------------------------------
    #toolbox.register("evaluate", fo.fitness)
    # fo.fitness2_1 e' il nuovo fitness con due elementi: impatto dulle deregolate e non deregolate
    #toolbox.register("evaluate", fo.fitness2_1)
    #toolbox.register("evaluate", getattr(fo, fitness))  ## NON FUNZIONA
    nomeFitness = getattr(fo, fitness)
    toolbox.register("evaluate", nomeFitness)
    #---------------------------------------------------------------------------

    toolbox.register("mate", cxTwoPointCopy)
    toolbox.register("mutate", tools.mutFlipBit, indpb=PROB_MUTAZIONE_SINGOLO_BIT)
    toolbox.register("select", tools.selTournament, tournsize=TOURNMENT_SIZE)

    pop = toolbox.population(n=POPOLAZIONE_INIZIALE)

    # Numpy equality function (operators.eq) between two arrays returns the
    # equality element wise, which raises an exception in the if similar()
    # check of the hall of fame. Using a different equality function like
    # numpy.array_equal or numpy.allclose solve this issue.
    hof = tools.HallOfFame(1, similar=np.array_equal)

    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    # tempo totale
    start_TOTALE = datetime.datetime.now().replace(microsecond=0)

    for i in range(NUM_RIPETIZIONI):
        # inizializza il cronometro
        start = datetime.datetime.now().replace(microsecond=0)
        # inizializza la popolazione
        pop = toolbox.population(n=POPOLAZIONE_INIZIALE)

        algorithms.eaSimple(pop,
                            toolbox,
                            cxpb=PROB_MATING,
                            mutpb=PROB_MUTAZIONE_INDIVIDUO,
                            ngen=NGEN,
                            stats=stats,
                            halloffame=hof)

        best_ind = tools.selBest(pop, 1)[0]
        # print "dim individuo", len(best_ind)
        print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values))
        print("\n")
        print("Iterazione " + str(i) + " di " + str(NUM_RIPETIZIONI))
        print("Scrittura su file\n")

        l_scelta = []
        for ii in range(len(best_ind)):
            if best_ind[ii] > 0:
                # memorizza il nome del miRNA
                l_scelta.append(l_miRNA_totale[ii])


        # cerco se ci sono proteine non coperte dai miRNA
        l_prot_non_coperte=protDeregNonCoperte(best_ind, l_p_deregolate)

        # SCRITTURA FILE #######################################################
        pre = "# "
        separatore = "-------------------------------------------------------\n"

        fout = open(OUTPUT_FILENAME, "a")

        fout.write("\n\n\n")
        fout.write(pre+"INIZIO RISULTATI" + separatore)

        tempoImpiegato = str( datetime.datetime.now().replace(microsecond=0) - start)
        stringa = pre+"fine processo: tempo impiegato "+ tempoImpiegato + "\n"
        fout.write(stringa)
        fout.write(pre + separatore)

        stringa = pre +"fitness migliore individuo " + str(best_ind.fitness.values) + "\n"
        fout.write(stringa)

        fout.write(pre + "ELENCO miRNA" + separatore)
        for item in l_scelta:
            stringa=pre + "> " + item +"\n"
            fout.write(stringa)



        fout.write(pre + "ELENCO PROTEINE NON COPERTE" + separatore)
        for item in l_prot_non_coperte:
            stringa=pre + "* " + item +"\n"
            fout.write(stringa)


        # CHIUSURA FILE DI OUTPUT
        fout.close()


    # Sccrive il tempo totale impiegato
    fout = open(OUTPUT_FILENAME, "a")
    tempoImpiegato_TOTALE = str( datetime.datetime.now().replace(microsecond=0) - start_TOTALE)
    stringa = "\n"+pre+"fine programma: tempo impiegato "+ tempoImpiegato_TOTALE + "\n"
    fout.write(stringa)
    fout.close()


    # manda una email per la fine dei calcoli
    # email risultati
    #---------------------------------------------------------------------------
    #subject = "Esperimenti "+ OUTPUT_FILENAME
    #body = "finiti tutti gli esperimenti: tempo impiegato "+tempoImpiegato
    #libPar.send_email(dm.user, dm.pwd, dm.recipient, subject, body)

    return pop, stats, hof, l_scelta

if __name__ == "__main__":

    pop, stats, hof, l_scelta=main(sys.argv[1])
