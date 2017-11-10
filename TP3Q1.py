import argparse

def buildKmer(k,input):
    kmerList=[]
    kmerPos=[]
    if k>0:
        # Parcourt la sequence
        for i in range(len(input) - k + 1):
            temp = ""
            # On ajoute les caracteres aux positions jusqu'a k
            for j in range(k):
                temp += input[i + j]
            # On ajoute chaque kmer individuellement
            kmerList.append(temp)
	    kmerPos.append(i)
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
        temp+="1"
    return temp

def findHSP(kmerList,seqDB, seed):
    compteurPos=0

    #Pour chaque sous-mot
    for kmer in kmerList:
	compteurPos = 0
	compteurSeed = 0
        #On loop sur la longueur de la sequence
        for compteurSeq in range(len(seqDB)):
            if seed[compteurSeed] == 1:
		    #Si le caractere est le meme dans le sous mot et la sequence
		    if kmer[compteurPos] == seqDB[compteurSeq]:
			#On continue d'incrementer dans le sous mot
			compteurPos += 1
		    #Sinon, on verifie le premier caractere pour repartir une sous boucle
		    elif kmer[0] == seqDB[compteurSeq]:
			charPos = 1
		    #Si rien ne match, on repart du premier caractere
		    else :
			charPos = 0
		    #On incremente la sequence totale
		    compteurSeq += 1
	    else:
		    #On a un zero, free pass
		    compteurSeq += 1
		    charPos += 1
	    compteurSeed+=1

def main():
    k = 5
    seqInput = "ACTGAAAATGAG"
    seqDB = "AGCATGACTGAAGTGAG"
    print(buildKmer(k, seqInput))
    seed = makeSeed(k)

if __name__ == "__main__":
    main()
