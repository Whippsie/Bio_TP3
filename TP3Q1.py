import argparse

seuil = 4


def buildKmer(k, inputUser):
    kmerList = []
    kmerPos = []
    if k > 0:
        # Parcourt la sequence
        for i in range(len(inputUser) - k + 1):
            temp = ""
            # On ajoute les caracteres aux positions jusqu'a k
            for j in range(k):
                temp += inputUser[i + j]
            # On ajoute chaque kmer individuellement
            kmerList.append(temp)
            kmerPos.append(i)
    return kmerList, kmerPos


def makeParser():
    parser = argparse.ArgumentParser(description='TP3_Bio')
    parser.add_argument('-i', help='Sequence ')
    parser.add_argument('-db', help='Base de donnees')
    parser.add_argument('-E', help='todo')
    parser.add_argument('-ss', help='todo')
    parser.add_argument('-seed', help='Graine 1 pour present, 0 pour no care')

    args = parser.parse_args()
    # print(args.accumulate(args.integers))


def makeSeed(k):
    temp = ""
    for i in range(k):
        temp += "1"
    return temp


def findHSP(kmerList, seqDB, seed):
    hspList = []
    hspPos = []
    # Pour chaque sous-mot
    for kmer in kmerList:
        compteurPos = 0
        # On loop sur la longueur de la sequence
        for compteurSeq in range(len(seqDB)):
            if seed[compteurPos] == "1":
                # Si le caractere est le meme dans le sous mot et la sequence
                if kmer[compteurPos] == seqDB[compteurSeq]:
                    if compteurPos == len(kmer) - 1:
                        # On a un mot!
                        hspList.append(kmer)
                        hspPos.append(compteurSeq)
                        compteurPos = 0
                        for i in range(len(kmer)):
                            if kmer[i] == seqDB[compteurSeq - (len(kmer) - 2)]:
                                compteurPos += 1
                                compteurSeq += 1
                            else:
                                compteurSeq += 1
                                compteurPos = 0
                    # On continue d'incrementer dans le sous mot
                    else:
                        compteurPos += 1
                        compteurSeq += 1
                # Sinon, on verifie le premier caractere pour repartir une sous boucle
                elif kmer[0] == seqDB[compteurSeq]:
                    compteurPos = 1
                # Si rien ne match, on repart du premier caractere
                else:
                    compteurPos = 0
                    # On incremente la sequence totale
                    compteurSeq += 1
            else:
                # On a un zero, free pass
                compteurSeq += 1
                compteurPos += 1
    print(hspList)
    print(hspPos)
    return (hspList, hspPos)


def extendGlouton(kmerList, hspList, kmerPos, seqDB, seq, hspPos):
    hspExtendedList = []
    for hsp in hspList:
        rightString = ""
        leftString = ""
        posDB = hspList.index(hsp)
        posDB = hspPos[posDB]
        maxScore = len(hsp) * 5
        currentScore = maxScore
        currentHSP = hsp
        while not bellowSeuil(maxScore, currentScore and not isEnd(currentHSP, seq, seqDB)):
            var = kmerList.index(hsp)
            pos = kmerPos[var]
            bothSidesString = currentHSP
            if pos > 0:
                bothSidesString = seq[pos - 1] + currentHSP
                tempChar = seq[pos - 1]
                tempCharDB = seqDB[posDB - 1]
                if tempChar == tempCharDB:
                    scoreLeft = 5 + currentScore
                    leftString = seq[pos - 1] + currentHSP
                else:
                    scoreLeft = currentScore - 4
                    leftString = currentHSP
            else:
                scoreLeft = -9000

            if posDB + len(currentHSP) < len(seqDB) and pos + len(currentHSP) < len(seq):
                bothSidesString += seq[pos + len(currentHSP)]
                tempChar = seq[pos + len(currentHSP)]
                tempCharDB = seqDB[posDB + len(currentHSP)]
                if tempChar == tempCharDB:
                    scoreRight = 5 + currentScore
                    rightString = currentHSP + seq[pos + len(currentHSP)]
                else:
                    scoreRight = currentScore - 4
                    rightString = currentHSP
            else:
                scoreRight = -9000

            scoreBoth = currentScore + (scoreLeft - currentScore) + (scoreRight - currentScore)
            currentScore = max(scoreLeft, scoreRight, scoreBoth)

            if currentScore == scoreBoth:
                currentHSP = bothSidesString
                pos -= 1
                posDB -= 1
            elif currentScore == scoreRight:
                currentHSP = rightString
            else:
                currentHSP = leftString
                pos -= 1
                posDB -= 1

            if currentScore > maxScore:
                maxScore = currentScore
        hspExtendedList.append(currentHSP)
    print(hspExtendedList)
    return hspExtendedList


def isEnd(currentHSP, seq, seqDB):
    if len(currentHSP) == len(seq) or len(currentHSP) == len(seqDB):
        return True
    return False


def bellowSeuil(maxScore, currentScore):
    return (maxScore - seuil) > currentScore


def main():
    k = 5
    seqInput = "ACTGAAAATGAG"
    seqDB = "AGCATGACTGAAGTGAG"
    (kmerList, kmerPos) = buildKmer(k, seqInput)
    print (kmerList)
    seed = makeSeed(k)
    (hspList, hspPos) = findHSP(kmerList, seqDB, seed)
    hspExtendedList = extendGlouton(kmerList, hspList, kmerPos, seqDB, seqInput, hspPos)


if __name__ == "__main__":
    main()
