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
            if temp not in kmerList:
                kmerList.append(temp)
                kmerPos.append(i)
    return kmerList, kmerPos


def makeParser():
    parser = argparse.ArgumentParser(description='TP3_Bio')
    parser.add_argument('-i', help='Sequence ')
    parser.add_argument('-db', help='Base de donnees')
    parser.add_argument('-E', help='E Value')
    parser.add_argument('-ss', help='SS')
    parser.add_argument('-seed', help='Graine 1 pour present, 0 pour no care')
    return parser.parse_args()


def makeSeed(k):
    temp = ""
    for i in range(k):
        temp += "1"
    return temp


def findHSP(kmerList, seqDB, seed):
    hspList = []
    hspPos = []
    infinite = 0
    # Pour chaque sous-mot
    for kmer in kmerList:
        compteurPos = 0
        # On loop sur la longueur de la sequence
        compteurSeq = 0
        while compteurSeq < len(seqDB):
            if seed[compteurPos] == "1":
                # Si le caractere est le meme dans le sous mot et la sequence
                if kmer[compteurPos] == seqDB[compteurSeq]:
                    if compteurPos == len(kmer) - 1:
                        # On a un mot!
                        hspList.append(kmer)
                        hspPos.append(compteurSeq - len(kmer) + 1)
                        if compteurSeq == (len(seqDB)-1):
                            compteurSeq += 1
                        else:
                            temp = (compteurSeq - compteurPos) + 1
                            if temp == (compteurSeq-len(kmer))+1:
                                compteurSeq = temp + 1
                            else:
                                compteurSeq = temp
                            compteurPos = 0
                        # On continue d'incrementer dans le sous mot
                    else:
                        compteurPos += 1
                        compteurSeq += 1
                # Sinon, on verifie le premier caractere pour repartir une sous boucle
                else:
                    if compteurPos > 0:
                        if compteurSeq == (len(seqDB)-1):
                            compteurSeq += 1
                        else:
                            if compteurPos > 1:
                                temp = (compteurSeq - compteurPos) + 1
                                if temp == (compteurSeq - len(kmer)) + 1:
                                    compteurSeq = temp + 1
                                else:
                                    compteurSeq = temp
                            compteurPos = 0
                    else:
                        compteurSeq += 1
            else:
                # On a un zero, free pass
                compteurPos += 1
                compteurSeq += 1
    print("hspList:", hspList)
    print("hspPos in original seq:", hspPos)
    return hspList, hspPos


def extendGlouton(kmerList, hspList, kmerPos, seqDB, seq, hspPos):
    hspExtendedList = []
    for hsp in hspList:
        rightString = ""
        leftString = ""
        #TODO: Probleme si plus qu'un HSP dans la sequence, nous donnera toujours le meme
        posDB = hspList.index(hsp)
        posDB = hspPos[posDB]
        maxScore = len(hsp) * 5
        currentScore = maxScore
        currentHSP = hsp
        pos = kmerList.index(hsp)
        pos = kmerPos[pos]

        while not bellowSeuil(maxScore, currentScore) and not isEnd(currentHSP, seq, seqDB):

            bothSidesString = currentHSP
            if pos > 0 and posDB > 0:
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

        hspExtendedList.append(Hsp(currentHSP, pos, pos + len(currentHSP) - 1, posDB, posDB + len(currentHSP) - 1))
    for hsp in hspExtendedList:
        print ("hspExtended:", hsp.hspString)
    return hspExtendedList


def isEnd(currentHSP, seq, seqDB):
    if len(currentHSP) == len(seq) or len(currentHSP) == len(seqDB):
        return True
    return False


def bellowSeuil(maxScore, currentScore):
    return (maxScore - seuil) >= currentScore

def merge(hspExtendedList, seqInput):
    i = 0
    while i < (len(hspExtendedList) - 1):
        for j in range((i+1), len(hspExtendedList)):
            hsp1 = hspExtendedList[i]
            hsp2 = hspExtendedList[j]
            if (hsp1.seqStart - hsp2.seqStart) == (hsp1.dbStart - hsp2.dbStart):
                newSeqStart = min([hsp1.seqStart, hsp2.seqStart])
                newSeqEnd = max([hsp1.seqEnd, hsp2.seqEnd])
                newDbStart = min([hsp1.dbStart, hsp2.dbStart])
                newDbEnd = max([hsp1.dbEnd, hsp2.dbEnd])
                newHspString = seqInput[newSeqStart:(newSeqEnd + 1)]
                hspExtendedList[i] = Hsp(newHspString, newSeqStart, newSeqEnd, newDbStart, newDbEnd)
                del hspExtendedList[j]
                break
            if j == len(hspExtendedList) - 1:
                i = i + 1
    for hsp in hspExtendedList:
        print ("hspMerged:", hsp.hspString)
    return hspExtendedList

class Hsp:
    def __init__(self, hspString, seqStart, seqEnd, dbStart, dbEnd):
        self.hspString = hspString
        self.seqStart = seqStart
        self.seqEnd = seqEnd
        self.dbStart = dbStart
        self.dbEnd = dbEnd
        
def calcBitScore(scoreBrut):
    lambd = 0.192
    k = 0.176
    return round(((lambd*scoreBrut)-numpy.log(k)/(numpy.log(2))), 6)


def calcEValue(lseq1, lseqdb, bitscore):
    return lseq1*lseqdb*(numpy.power(2, -bitscore))

def main():
    k = 4
    seqInput = "ACGCGCGAAGAGCG"
    seqDB = "TACGCGCGAAGCG"
    kmerList, kmerPos = buildKmer(k, seqInput)
    print ("kmerList:", kmerList)
    seed = makeSeed(k)
    hspList, hspPos = findHSP(kmerList, seqDB, seed)
    hspExtendedList = extendGlouton(kmerList, hspList, kmerPos, seqDB, seqInput, hspPos)
    hspMergedList = merge(hspExtendedList, seqInput)
    args = makeParser()
    #NOTE : Args all treated as string
    #if args.i:
    #.format(args.square, answer)
    # print(args.accumulate(args.integers))

if __name__ == "__main__":
    main()
