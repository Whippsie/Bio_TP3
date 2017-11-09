import argparse

def buildKmer(k,input):
    kmerList=[]
    if k>0:
        # Parcourt la sequence
        for i in range(len(input) - k + 1):
            temp = ""
            # On ajoute les caracteres aux positions jusqu'a k
            for j in range(k):
                temp += input[i + j]
            # On ajoute chaque kmer individuellement
            kmerList.append(temp)
    return kmerList

def makeParser():
    parser = argparse.ArgumentParser(description='TP3_Bio')
    parser.add_argument('-i', help='Sequence ')
    parser.add_argument('-db', help='Base de donnees')
    parser.add_argument('-E', help='todo')
    parser.add_argument('-ss', help='todo')
    parser.add_argument('-seed', help='Graine 1 pour present, 0 pour no care')

    args = parser.parse_args()
    #print(args.accumulate(args.integers))

def makeSeed(k):
    temp=""
    for i in range(k):
        temp+=1
    return temp

def findHSP(kmerList,seqDB):
    compteurSeq=0
    compteurKmer=0

    #TODO: Rajoute le loop dans la sequence
    #Pour chaque sous-mot
    for kmer in kmerList:
        #On loop sur la longueur du sous-mot
        for char in range(len(kmer)):
            #TODO:Rajouter la prise en charge du 1 et du 0 dans le seed
            #Si le caractere est le meme dans le sous mot et la sequence
            if char == seqDB[compteurSeq]:
                #On continue d'incrementer dans le sous mot
                char += 1
            else:
                #On ressaiera a partir du debut au prochain caractere
                char = 0
            #On incremente la sequence totale
            compteurSeq += 1

def main():
    k = 5
    seqInput = "ACTGAAAATGAG"
    seqDB = "AGCATGACTGAAGTGAG"
    print(buildKmer(k, seqInput))
    seed = makeSeed(k)

if __name__ == "__main__":
    main()