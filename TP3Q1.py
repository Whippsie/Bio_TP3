import argparse

seuil = 0

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
		#TODO: Changer pour pos commencement
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
	hspList=[]
	hspPos=[]
	compteurPos=0

    #Pour chaque sous-mot
	for kmer in kmerList:
		compteurPos = 0
		#On loop sur la longueur de la sequence
		for compteurSeq in range(len(seqDB)):
			if seed[compteurPos] == "1":
		    #Si le caractere est le meme dans le sous mot et la sequence
				if kmer[compteurPos] == seqDB[compteurSeq]:
					if compteurPos == len(kmer) - 1:
						#On a un mot!
						hspList.append(kmer)
						hspPos.append(compteurSeq)
						compteurPos = 0
						for i in range(len(kmer)):
							if kmer[i] == seqDB[compteurSeq-(len(kmer)-2)]:						
								compteurPos += 1
								compteurSeq += 1
							else:
								compteurSeq += 1
								compteurPos = 0
					#On continue d'incrementer dans le sous mot
					else:
						compteurPos += 1
						compteurSeq += 1
				#Sinon, on verifie le premier caractere pour repartir une sous boucle
				elif kmer[0] == seqDB[compteurSeq]:
					compteurPos = 1
					#Si rien ne match, on repart du premier caractere
				else :
					compteurPos = 0
					#On incremente la sequence totale
					compteurSeq += 1
			else:
				#On a un zero, free pass
				compteurSeq += 1
				compteurPos += 1
	print(hspList)
	print(hspPos)


def extendGlouton(kmerList, hspList, kmerPos, seqDB, seq, hspPos):
	for hsp in hspList:
		posDB = hspList.index[hsp]
		posDB = hspPos[posDB]
		maxScore = len(hsp)*5
		currentScore = maxScore
		currentHSP = hsp		
		while(!bellowSeuil(maxScore, currentScore) and !isEnd()):
			var = kmerList.index[hsp]
			pos = kmerPos[var]
			scoreLeft = 0
			scoreRight = 0
			if pos > 0:
				tempChar = seq[pos-1]
				tempCharDB = seqDB[posDB-1]
				if tempChar == tempCharDB:
					scoreLeft = 5 + currentScore
				else:
					scoreLeft = currentScore - 4
			if pos+len(currentHSP) < len(seqDB) or pos+len(currentHSP) < len(seq):
				tempChar = seq[pos+1]
				tempCharDB = seqDB[posDB+1]
				if tempChar == tempCharDB:
					scoreRight = 5 + currentScore
				else:
					scoreRight = currentScore - 4
			scoreBoth = currentScore + (scoreLeft - currentScore) + (scoreRight - currentScore)
			currentScore = max(scoreLeft, scoreRight, scoreBoth)
			if currentScore == scoreBoth:
				currentHSP = seq[pos-1] + currentHSP + seq[pos+1]
			elif currentScore == scoreRight:
				currentHSP = currentHSP + seq[pos+1]
			else:
				currentHSP = seq[pos-1] + currentHSP
			if currentScore > maxScore:
				maxScore = currentScore

def getMax(sl,sr,sb,posDB,pos,currHSP):
	if sl>sr:
		if sl>sb:
			
			return (sl, posDB-1, pos-1)
		else:
			return (sb, 

def isEnd():
	return False

def bellowSeuil(maxScore, val):
	return (maxScore-seuil) <= val


def main():
	k = 5
	seqInput = "ACTGAAAATGAG"
	seqDB = "AGCATGACTGAAGTGAG"
	kmerList = buildKmer(k, seqInput)
	print (kmerList)
	seed = makeSeed(k)
	findHSP(kmerList,seqDB,seed)

if __name__ == "__main__":
	main()
