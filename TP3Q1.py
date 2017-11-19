#!/usr/bin/env python
# encoding: utf-8

import argparse
import math

seuil = 12
seuilCutoff = 0.1
match, mismatch, indel = 1, -1, -1


def fastaSequences(filename):
    sequences = []
    with open(filename, 'r') as f:
        for line in f:
            sequences.append(line.strip("\n"))
    return sequences


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
    hspPosDB = []
    hspPosOriginale = []
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
                        hspPosDB.append(compteurSeq - len(kmer) + 1)
                        if compteurSeq == (len(seqDB) - 1):
                            compteurSeq += 1
                        else:
                            temp = (compteurSeq - compteurPos) + 1
                            if temp == (compteurSeq - len(kmer)) + 1:
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
                        if compteurSeq == (len(seqDB) - 1):
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
    # print("hspList:", hspList)
    # print("hspPos in DB seq:", hspPosDB)
    return hspList, hspPosDB


def extendGlouton(kmerList, hspList, kmerPos, seqDB, seq, hspPos):
    hspExtendedList = []
    hspScoreList = []
    i = 0
    for hsp in hspList:
        rightString = ""
        leftString = ""
        # i is the position of the current hsp, we don't use index to protect the twin hsp who wouldn't be caught
        posDB = hspPos[i]
        maxScore = len(hsp) * 5
        currentScore = maxScore
        currentHSP = hsp
        lastScore = maxScore
        lastHSP = hsp
        #TODO: PAS bon si 2 fois le même hsp dans les kmers
        pos = kmerList.index(hsp)
        pos = kmerPos[pos]
        while not bellowSeuil(maxScore, currentScore) and not isEnd(currentHSP, seq, seqDB):
            bothSidesString = currentHSP
            # print("pos:", pos," posDB:", posDB)
            if pos > 0 and posDB > 0:
                bothSidesString = seq[pos - 1] + currentHSP
                tempChar = seq[pos - 1]
                tempCharDB = seqDB[posDB - 1]
                if tempChar == tempCharDB:
                    scoreLeft = 5 + currentScore
                    leftString = seq[pos - 1] + currentHSP
                else:
                    scoreLeft = currentScore - 4
                    leftString = seq[pos - 1] + currentHSP
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
                    rightString = currentHSP + seq[pos + len(currentHSP)]
            else:
                scoreRight = -9000

            scoreBoth = currentScore + (scoreLeft - currentScore) + (scoreRight - currentScore)
            currentScore = max(scoreLeft, scoreRight, scoreBoth)

            if currentScore == scoreBoth or currentScore == -9000:
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

            if not bellowSeuil(maxScore, currentScore):
                lastHSP = currentHSP
                lastScore = currentScore
            else:
                break

        hspExtendedList.append(Hsp(lastHSP, pos, pos + len(currentHSP) - 1, posDB, posDB + len(currentHSP) - 1, lastScore))
        hspScoreList.append(lastScore)
        i += 1
    #for hsp in hspExtendedList:
        #print ("hspExtended:", hsp.hspString)
    return hspExtendedList, hspScoreList


def isEnd(currentHSP, seq, seqDB):
    if len(currentHSP) == len(seq) or len(currentHSP) == len(seqDB):
        return True
    return False


def bellowSeuil(maxScore, currentScore):
    return (maxScore - currentScore) >= seuil


def merge(hspExtendedList, seqInput, seqDB):
    i = 0
    while i < (len(hspExtendedList) - 1):
        for j in range((i + 1), len(hspExtendedList)):
            hsp1 = hspExtendedList[i]
            hsp2 = hspExtendedList[j]
            if mergeable(hsp1, hsp2):
                newSeqStart = min([hsp1.seqStart, hsp2.seqStart])
                newSeqEnd = max([hsp1.seqEnd, hsp2.seqEnd])
                newDbStart = min([hsp1.dbStart, hsp2.dbStart])
                newDbEnd = max([hsp1.dbEnd, hsp2.dbEnd])
                newHspString = seqInput[newSeqStart:(newSeqEnd + 1)]
                newScore = calcHspBruteScore(newHspString, seqDB, newDbStart)
                hspExtendedList[i] = Hsp(newHspString, newSeqStart, newSeqEnd, newDbStart, newDbEnd, newScore)
                del hspExtendedList[j]
                break
            if j == len(hspExtendedList) - 1:
                i = i + 1
    return hspExtendedList

def mergeable(hsp1, hsp2):
    if (hsp1.dbStart <= hsp2.dbStart and hsp1.dbEnd >= hsp2.dbStart) or (hsp1.dbStart <= hsp2.dbEnd and hsp1.dbEnd >= hsp2.dbEnd):
        #print("almost there", hsp1.hspString, "|", hsp2.hspString)
        if (hsp1.seqStart - hsp2.seqStart) == (hsp1.dbStart - hsp2.dbStart):
            #print("mergeable",hsp1.hspString,"|",hsp2.hspString)
            return True
    return False

def calcHspBruteScore(hspString, seqDB, dbStart):
    newScore = 0
    j = dbStart
    for i in range(0, len(hspString)):
        if hspString[i] == seqDB[j]:
            newScore += 5
        else:
            newScore += -4
        j += 1
    return newScore


class Hsp:
    def __init__(self, hspString, seqStart, seqEnd, dbStart, dbEnd, score):
        self.hspString = hspString
        self.seqStart = seqStart
        self.seqEnd = seqEnd
        self.dbStart = dbStart
        self.dbEnd = dbEnd
        self.score = score

class HspForFilter:
    def __init__(self, hsp, bitscore, eValue, score):
        self.hsp = hsp
        self.bitscore = bitscore
        self.eValue = eValue
        self.score = score

def filterHSP(hspList, l1, ldb):
    hspWithScores = []
    for hsp in hspList:
        bitscore = calcBitScore(hsp.score)
        eValue = calcEValue(l1, ldb, bitscore)
        hspWithScores.append(HspForFilter(hsp, bitscore, eValue, hsp.score))
    bestBitScore = -9000
    selectedHsp = None

    for hsp in hspWithScores:
        if hsp.eValue <= seuilCutoff and hsp.bitscore > bestBitScore:
            bestBitScore = hsp.bitscore
            selectedHsp = hsp

    return selectedHsp


def calcBitScore(scoreBrut):
    lambd = 0.192
    k = 0.176
    return round((((lambd * scoreBrut) - math.log(k)) / (math.log(2))), 0)


def calcEValue(lseq1, lseqdb, bitscore):
    return lseq1 * lseqdb * (math.pow(2, -bitscore))


def getLengthSeqDB(seqSearchDB):
    lgt = 0
    for seq in seqSearchDB:
        if '>' not in seq:
            lgt += len(seq)
    return lgt

#Code Alignement

class HspAlignment:
    def __init__(self, seqInput, seqDB, score, inputStart, inputEnd, dbStart, dbEnd):
        self.seqInput = seqInput
        self.seqDB = seqDB
        self.score = score
        self.inputStart = inputStart
        self.inputEnd = inputEnd
        self.dbStart = dbStart
        self.dbEnd = dbEnd

def alignment(seq1, seq2):
    #print(seq1)
    #print(seq2)
    matrix, directions = createMatrix(seq1, seq2)
    #printMatrix(matrix)
    #printMatrix(directions)
    solvedMatrix, directionMatrix, maxPos, maxScore = solveMatrix(matrix, directions, seq1, seq2)
    #printMatrix(solvedMatrix)
    #printMatrix(directionMatrix)
    hspAlignment = traceback(directionMatrix, seq1, seq2, maxPos, maxScore)
    return hspAlignment


def createMatrix(seq1, seq2):
    scoreMatrix = []
    directionMatrix = []
    rows = len(seq1) + 1
    columns = len(seq2) + 1

    for row in range(rows):
        rowArrayScore = []
        rowArrayDir = []
        for column in range(columns):
            rowArrayScore.append(0)
            rowArrayDir.append(0)
        scoreMatrix.append(rowArrayScore)
        directionMatrix.append(rowArrayDir)

    return scoreMatrix, directionMatrix

def solveMatrix (matrix, directionMatrix, seq1, seq2):
    rows = len(seq1) + 1
    columns = len(seq2) + 1
    maxScore = 0
    maxPos = None

    for i in range(1, rows):
        for j in range(1, columns):
            score, dir = evalScore(matrix, i, j, seq1, seq2)
            matrix[i][j] = score
            directionMatrix[i][j] = dir
            if score > maxScore:
                maxScore = score
                maxPos = (i, j)

    return matrix, directionMatrix, maxPos, maxScore

def evalScore (matrix, i, j, seq1, seq2):
    if seq1[i - 1] == seq2[j - 1]:
        diag = matrix[i - 1][j - 1] + match
    else:
        diag = matrix[i - 1][j - 1] + mismatch

    up = matrix[i-1][j] + indel
    left = matrix[i][j-1] + indel

    score = max(0, diag, up, left)

    #Directions:
    #0 = fin
    #1 = diagonale
    #2 = haut
    #3 = gauche
    dir = 0
    if diag == score:
        dir = 1
    elif up == score:
        dir = 2
    elif left == score:
        dir = 3
    return score, dir

def printMatrix(matrix):
    print('\n'.join([''.join(['{:4}'.format(item) for item in row])
                     for row in matrix]))

def traceback (directions, seq1, seq2, start, score):
    i, j = start
    inputEnd = i - 1
    dbEnd = j - 1
    trace = directions[i][j]
    seq1Aligned = []
    seq2Aligned = []

    while trace != 0:
        if trace == 1:
            seq1Aligned.append(seq1[i-1])
            seq2Aligned.append(seq2[j-1])
            i -= 1
            j -= 1
        elif trace == 2:
            seq1Aligned.append(seq1[i-1])
            seq2Aligned.append('-')
            i -= 1
        else:
            seq2Aligned.append(seq2[j-1])
            seq1Aligned.append('-')
            j -= 1

        trace = directions[i][j]

    if (i-1) >= 0 and (j-1) >= 0:
        seq1Aligned.append(seq1[i - 1])
        inputStart = i - 1
        seq2Aligned.append(seq2[j - 1])
        dbStart = j - 1
    else:
        inputStart = i
        dbStart = j

    seq1Aligned = ''.join(reversed(seq1Aligned))
    seq2Aligned = ''.join(reversed(seq2Aligned))
    print("Meilleur alignement")
    print ("Input:",seq1Aligned)
    print ("DB:   ",seq2Aligned)
    print("suppose etre")
    print("GAAAATCCTCGTGTCACCAGTTCAAATCTGGTTCCTGGCA")
    print("GAAAATCCTTGTGTCAGTGGTTCAAATCCACTTTCAGGCA")
    hspAlignment = HspAlignment(seq1Aligned, seq2Aligned, score, inputStart, inputEnd, dbStart, dbEnd)

    return hspAlignment

def printAlignment (hspAlignment):
    dbStart = fixStartIndices(str(hspAlignment.dbStart))
    dbEnd = fixEndIndices(str(hspAlignment.dbEnd))
    inputStart = fixStartIndices(str(hspAlignment.inputStart))
    inputEnd = fixEndIndices(str(hspAlignment.inputEnd))

    print ("Alignement:")
    print (inputStart + hspAlignment.seqInput + inputEnd)
    print (dbStart + hspAlignment.seqDB + dbEnd)
    print ("---------------------------------------------")
    #print ("Score: ", hspAlignment.score)

def fixStartIndices (str):
    while len(str) < 3:
        str += " "
    return str

def fixEndIndices (str):
    while len(str) < 3:
        str = " " + str
    return str
#Code Alignement

def main():
    # seed = makeSeed(k)
    # POUR TESTER UNIQUEMENT
    seed = "110011"
    #POUR TESTER UNIQUEMENT
    k = len(seed)
    seqSearch = fastaSequences("unknown.fasta")
    #seqInput = seqSearch[0]

    # tests
    # seqInput = "ACGCGCGAAGAGCG"
    # seqDB = "TACGCGCGAAGCG"
    # seqDB = "TACGCGTGAAGCG"
    # seqDB = "TACGCGCGAAACG"
    # seqDB = "TACGCGTGAAACG"

    seqSearchDB = fastaSequences("tRNA2s.fasta")
    #seqDB = seqSearchDB[0]

    #for seqInput in seqSearch:
    seqInput = "GAAAATCCTCGTGTCACCAGTTCAAATCTGGTTCCTGGCA"
    for seqDB in seqSearchDB:
        if '>' in seqDB:
            print(seqDB)
            temp = seqDB
            continue
        # MARCHE PAS POUR L'INSTANT
        # alignLocal(seqInput, seqDB)
        kmerList, kmerPos = buildKmer(k, seqInput)


        hspList, hspPos = findHSP(kmerList, seqDB, seed)
       # print("list:",hspList)
        hspExtendedList, hspScoreList = extendGlouton(kmerList, hspList, kmerPos, seqDB, seqInput, hspPos)
        hspMergedList = merge(hspExtendedList, seqInput, seqDB)
        #for hsp in hspMergedList:
        #    print ("hsp merged",hsp.hspString)
        # PAS FINI
        dbSeqLength = getLengthSeqDB(seqSearchDB)
        selectedHsp = filterHSP(hspMergedList, len(seqInput), dbSeqLength)
        if selectedHsp is not None:
            bitscore = selectedHsp.bitscore
            eValue = selectedHsp.eValue
            #hspAlign = alignment(selectedHsp.hsp.hspString, seqDB)
            hspAlign = alignment(seqInput, seqDB)
            print (temp , " Score: ", hspAlign.score, " Ident: TODO")
            print("seq input : ", seqInput)
            print("seq DB : ", seqDB)
            print ("Selected hsp: ", selectedHsp.hsp.hspString)
            print (" ")
            print ("# Best HSP:")
            print ("Id:",temp," Score brut:", selectedHsp.score, " Bitscore:", bitscore," Evalue: ", eValue)
            printAlignment(hspAlign)
            print("Suppose etre")
            print("GAAAATCCTCGTGTCACCAGTTCAAATC")
            print("GAAAATCCTTGTGTCAGTGGTTCAAATC")
        args = makeParser()

        # NOTE : Args all treated as string
        # if args.i:
        # .format(args.square, answer)
        # print(args.accumulate(args.integers))

if __name__ == "__main__":
    main()
