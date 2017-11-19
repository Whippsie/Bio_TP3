#!/usr/bin/env python
# encoding: utf-8

import argparse
import math

seuil = 4
match, mismatch, indel = 1, -1, -1


def fastaSequences(filename):
    sequences = []
    with open(filename, 'r') as f:
        for line in f:
            if '>' not in line[0]:
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
    print("hspList:", hspList)
    print("hspPos in DB seq:", hspPosDB)
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

        hspExtendedList.append(Hsp(currentHSP, pos, pos + len(currentHSP) - 1, posDB, posDB + len(currentHSP) - 1))
        hspScoreList.append(maxScore)
        i += 1
    for hsp in hspExtendedList:
        print ("hspExtended:", hsp.hspString)
    return hspExtendedList, hspScoreList


def isEnd(currentHSP, seq, seqDB):
    if len(currentHSP) == len(seq) or len(currentHSP) == len(seqDB):
        return True
    return False


def bellowSeuil(maxScore, currentScore):
    return (maxScore - seuil) >= currentScore


def merge(hspExtendedList, seqInput):
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
                hspExtendedList[i] = Hsp(newHspString, newSeqStart, newSeqEnd, newDbStart, newDbEnd)
                del hspExtendedList[j]
                break
            if j == len(hspExtendedList) - 1:
                i = i + 1
    for hsp in hspExtendedList:
        print ("hspMerged:", hsp.hspString)
    return hspExtendedList

def mergeable(hsp1, hsp2):
    if (hsp1.dbStart <= hsp2.dbStart and hsp1.dbEnd >= hsp2.dbStart) or (hsp1.dbStart <= hsp2.dbEnd and hsp1.dbEnd >= hsp2.dbEnd):
        if (hsp1.seqStart - hsp2.seqStart) == (hsp1.dbStart - hsp2.dbStart):
            return True
    return False

def filterMergedList(hspMergedList, seqDB):
    hspMerged = ""
    maxScore = -9000
    for hsp in hspMergedList:
        tempHspAlignment = alignment(hsp.hspString, seqDB)
        if tempHspAlignment.score > maxScore:
            maxScore = tempHspAlignment.score
            hspMerged = tempHspAlignment.seqInput
    print("Chosen HSP: ", hspMerged)
    print("Alignment Score: ", maxScore)
    return hspMerged, maxScore

class Hsp:
    def __init__(self, hspString, seqStart, seqEnd, dbStart, dbEnd):
        self.hspString = hspString
        self.seqStart = seqStart
        self.seqEnd = seqEnd
        self.dbStart = dbStart
        self.dbEnd = dbEnd


def filterHSP(hspList, hspScoreList, l1, ldb):
    i = 0
    for hsp in hspList:
        bitscore = calcBitScore(hspScoreList[i])
        calcEValue(l1, ldb, bitscore)
        i += 1


def calcBitScore(scoreBrut):
    lambd = 0.192
    k = 0.176
    return round(((lambd * scoreBrut) - math.log(k) / (math.log(2))), 6)


def calcEValue(lseq1, lseqdb, bitscore):
    return lseq1 * lseqdb * (math.pow(2, -bitscore))


def getLengthSeqDB(seqSearchDB):
    lgt = 0
    for seq in seqSearchDB:
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
    matrix, directions = createMatrix(seq1, seq2)
    solvedMatrix, maxPos, maxScore = solveMatrix(matrix, directions, seq1, seq2)
    hspAlignment = traceback(directions, seq1, seq2, maxPos, maxScore)
    printAlignment(hspAlignment)
    print (solvedMatrix)
    return hspAlignment


def createMatrix(seq1, seq2):
    scoreMatrix = []
    directionMatrix = []
    rows = len(seq1) + 1
    columns = len(seq2) + 1

    for row in range(rows):
        rowArray = []
        for column in range(columns):
            rowArray.append(0)
        scoreMatrix.append(rowArray)
        directionMatrix.append(rowArray)

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
            print ("(i, j, score, dir)", i, j, score, dir)
            print ("=======")
            directionMatrix[i][j] = dir
            if score > maxScore:
                maxScore = score
                maxPos = (i, j)

    print (matrix)

    return matrix, maxPos, maxScore

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

    print ("(i, j, score, dir)", i, j, score, dir)
    return score, dir

def traceback (directions, seq1, seq2, start, score):
    i, j = start
    inputEnd = i - 1
    dbEnd = j - 1
    move = directions[i][j]
    seq1Aligned = []
    seq2Aligned = []

    while move != 0:
        if move == 1:
            seq1Aligned.append(seq1[i-1])
            seq2Aligned.append(seq2[j-1])
            i -= 1
            j -= 1
        elif move == 2:
            seq1Aligned.append(seq1[i-1])
            seq2Aligned.append('-')
            i -= 1
        else:
            seq2Aligned.append(seq2[j-1])
            seq1Aligned.append('-')
            j -= 1

        move = directions[i][j]

    if (i-1) >= 0:
        seq1Aligned.append(seq1[i - 1])
        inputStart = i - 1
    else:
        inputStart = i

    if (j-1) >= 0:
        seq2Aligned.append(seq2[j - 1])
        dbStart = j - 1
    else:
        dbStart = j

    seq1Aligned = ''.join(reversed(seq1Aligned))
    seq2Aligned = ''.join(reversed(seq2Aligned))

    hspAlignment = HspAlignment(seq1Aligned, seq2Aligned, score, inputStart, inputEnd, dbStart, dbEnd)

    return hspAlignment

def printAlignment (hspAlignment):
    dbStart = fixStartIndices(str(hspAlignment.dbStart))
    dbEnd = fixEndIndices(str(hspAlignment.dbEnd))
    inputStart = fixStartIndices(str(hspAlignment.inputStart))
    inputEnd = fixEndIndices(str(hspAlignment.inputEnd))

    print ("Alignement:")
    print (dbStart + hspAlignment.seqDB + dbEnd)
    print (inputStart + hspAlignment.seqInput + inputEnd)
    print ("Score: ", hspAlignment.score)

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
    print(" ")
    print(" ")
    print("DELIMITER====================================================================")
    k = 4
    seqSearch = fastaSequences("unknown.fasta")
    seqInput = seqSearch[0]

    # tests
    # seqInput = "ACGCGCGAAGAGCG"
    # seqDB = "TACGCGCGAAGCG"
    # seqDB = "TACGCGTGAAGCG"
    # seqDB = "TACGCGCGAAACG"
    # seqDB = "TACGCGTGAAACG"

    seqSearchDB = fastaSequences("tRNAs.fasta")
    seqDB = seqSearchDB[0]

    print("seq input : ", seqInput)
    print("seq DB : ", seqDB)

    # MARCHE PAS POUR L'INSTANT
    # alignLocal(seqInput, seqDB)
    kmerList, kmerPos = buildKmer(k, seqInput)
    seed = makeSeed(k)

    hspList, hspPos = findHSP(kmerList, seqDB, seed)
    hspExtendedList, hspScoreList = extendGlouton(kmerList, hspList, kmerPos, seqDB, seqInput, hspPos)
    hspMergedList = merge(hspExtendedList, seqInput)
    hspMerged, hspMergedScore = filterMergedList(hspMergedList, seqDB)

    # PAS FINI
    dbSeqLength = getLengthSeqDB(seqSearchDB)
    filterHSP(hspList, hspScoreList, len(seqInput), dbSeqLength)
    args = makeParser()

    # NOTE : Args all treated as string
    # if args.i:
    # .format(args.square, answer)
    # print(args.accumulate(args.integers))


if __name__ == "__main__":
    main()
