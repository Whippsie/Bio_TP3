#!/usr/bin/env python
# encoding: utf-8

import argparse
import math
import time

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
                kmerList.append(Kmer(temp, i))
    return kmerList


class Kmer:
    def __init__(self, kmerString, seqStart):
        self.kmerString = kmerString
        self.seqStart = seqStart
        self.seqEnd = seqStart + len(kmerString)


def findHSP(kmerList, seqDB, seed):
    hspList = []
    hspPosDB = []
    # Pour chaque sous-mot
    for kmerObj in kmerList:
        kmer = kmerObj.kmerString
        temp = ""
        for i in range(len(seqDB) - len(kmer)):
            j = 0
            while (j < len(kmer)) and (seed[j] == "0" or seqDB[i + j] == kmer[j]):
                j += 1
            temp += seqDB[i]
            if j == len(kmer):
                hspList.append(kmerObj)
                hspPosDB.append(i)
    return hspList, hspPosDB


def extendGlouton(hspList, seqDB, seq, hspPos, seuil):
    hspExtendedList = []
    i = 0
    for hspObj in hspList:
        hsp = hspObj.kmerString
        rightString = ""
        leftString = ""
        # i is the position of the current hsp, we don't use index to protect the twin hsp who wouldn't be caught
        posDB = hspPos[i]
        maxScore = len(hsp) * 5
        currentScore = maxScore
        currentHSP = hsp
        lastScore = maxScore
        lastHSP = hsp
        pos = hspObj.seqStart
        while not bellowSeuil(maxScore, currentScore, seuil) and not isEnd(currentHSP, seq, seqDB):
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
            elif currentScore == scoreLeft:
                currentHSP = leftString
                pos -= 1
                posDB -= 1
            else:
                currentHSP = rightString

            if currentScore > maxScore:
                maxScore = currentScore

            if not bellowSeuil(maxScore, currentScore, seuil):
                lastHSP = currentHSP
                lastScore = currentScore
            else:
                break
        temp = Hsp(lastHSP, pos, pos + len(currentHSP) - 1, posDB, posDB + len(currentHSP) - 1, lastScore)
        hspExtendedList.append(temp)
        i += 1
    return hspExtendedList


def isEnd(currentHSP, seq, seqDB):
    if len(currentHSP) == len(seq) or len(currentHSP) == len(seqDB):
        return True
    return False


def bellowSeuil(maxScore, currentScore, seuil):
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
    if (hsp1.dbStart <= hsp2.dbStart and hsp1.dbEnd >= hsp2.dbStart) or (
            hsp1.dbStart <= hsp2.dbEnd and hsp1.dbEnd >= hsp2.dbEnd):
        if (hsp1.seqStart - hsp2.seqStart) == (hsp1.dbStart - hsp2.dbStart):
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


def filterHSP(hspList, l1, ldb, seuilCutOff):
    hspWithScores = []
    for hsp in hspList:
        bitscore = calcBitScore(hsp.score)
        eValue = calcEValue(l1, ldb, bitscore)
        hspWithScores.append(HspForFilter(hsp, bitscore, eValue, hsp.score))
    bestBitScore = -9000
    selectedHsp = None

    for hsp in hspWithScores:
        if hsp.eValue <= seuilCutOff and hsp.bitscore > bestBitScore:
            bestBitScore = hsp.bitscore
            selectedHsp = hsp

    return selectedHsp


def calcBitScore(scoreBrut):
    lambd = 0.192
    k = 0.176
    return round((((lambd * scoreBrut) - math.log(k)) / (math.log(2))), 0)


def calcEValue(lseq1, lseqdb, bitscore):
    try:
        ans = lseq1 * lseqdb * (math.pow(2, -bitscore))
    except OverflowError:
        ans = float('inf')
    return ans


def getLengthSeqDB(seqSearchDB):
    lgt = 0
    for seq in seqSearchDB:
        if '>' not in seq:
            lgt += len(seq)
    return lgt


# Code Alignement
class HspAlignment:
    def __init__(self, seqInput, seqDB, score, inputStart, inputEnd, dbStart, dbEnd):
        self.seqInput = seqInput
        self.seqDB = seqDB
        self.score = score
        self.inputStart = inputStart
        self.inputEnd = inputEnd
        self.dbStart = dbStart
        self.dbEnd = dbEnd
        self.ident = ident(seqInput, seqDB)


def ident(seq1, seq2):
    sim = float(0)
    for i in range(0, len(seq1)):
        if seq1[i] == seq2[i]:
            sim += 1
    return sim / len(seq1)


def alignment(seq1, seq2):
    matrix, directions = createMatrix(seq1, seq2)
    solvedMatrix, directionMatrix, maxPos, maxScore = solveMatrix(matrix, directions, seq1, seq2)
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


def solveMatrix(matrix, directionMatrix, seq1, seq2):
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


def evalScore(matrix, i, j, seq1, seq2):
    if seq1[i - 1] == seq2[j - 1]:
        diag = matrix[i - 1][j - 1] + match
    else:
        diag = matrix[i - 1][j - 1] + mismatch

    up = matrix[i - 1][j] + indel
    left = matrix[i][j - 1] + indel

    score = max(0, diag, up, left)

    # Directions:
    # 0 = fin
    # 1 = diagonale
    # 2 = haut
    # 3 = gauche
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


def traceback(directions, seq1, seq2, start, score):
    i, j = start
    inputEnd = i - 1
    dbEnd = j - 1
    trace = directions[i][j]
    seq1Aligned = []
    seq2Aligned = []

    while trace != 0:
        if trace == 1:
            seq1Aligned.append(seq1[i - 1])
            seq2Aligned.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif trace == 2:
            seq1Aligned.append(seq1[i - 1])
            seq2Aligned.append('-')
            i -= 1
        else:
            seq2Aligned.append(seq2[j - 1])
            seq1Aligned.append('-')
            j -= 1

        trace = directions[i][j]

    if (i - 1) >= 0 and (j - 1) >= 0:
        seq1Aligned.append(seq1[i - 1])
        inputStart = i - 1
        seq2Aligned.append(seq2[j - 1])
        dbStart = j - 1
    else:
        inputStart = i
        dbStart = j

    seq1Aligned = ''.join(reversed(seq1Aligned))
    seq2Aligned = ''.join(reversed(seq2Aligned))
    hspAlignment = HspAlignment(seq1Aligned, seq2Aligned, score, inputStart, inputEnd, dbStart, dbEnd)

    return hspAlignment


def printSmithWaterman(hspAlignment):
    print (hspAlignment.seqDB)
    print (hspAlignment.seqInput)
    print (" ")


# Code found on Stackoverflow and modified to be of use here
def quicksort(array):
    less = []
    equal = []
    greater = []

    if len(array) > 1:
        pivot = array[0].selectedHsp.score
        for x in array:
            if x.selectedHsp.score < pivot:
                less.append(x)
            if x.selectedHsp.score == pivot:
                equal.append(x)
            if x.selectedHsp.score > pivot:
                greater.append(x)

        return quicksort(greater) + equal + quicksort(less)
    else:
        return array


class Result:
    def __init__(self, selectedHsp, seqDB, description):
        self.selectedHsp = selectedHsp
        self.seqDB = seqDB
        self.description = description


def printAlignment(result, selectedHsp, bitscore, eValue, dbSeqLength):
    seqDB = result.seqDB
    hspAlignment = selectedHsp.hsp
    seqHSP = hspAlignment.hspString
    posDB = hspAlignment.dbStart
    i=0
    if seqDB[posDB+1] != seqHSP[i] and seqDB[posDB]!=seqHSP[i]:
        if seqDB[posDB+1] == seqHSP[i+1] :
            # On doit shifter de 1 le HSP
            temp = ""
            for j in range(len(seqHSP)):
                if j != 0:
                    temp += seqHSP[j]
            seqHSP = temp
            hspAlignment.dbStart+=1
            hspAlignment.dbEnd+=1
            hspAlignment.seqStart+=2
            hspAlignment.seqEnd+=2
    elif seqDB[posDB] != seqHSP[i]:
        #On coupe les deux, mismatch
        selectedHsp.score += 4
        bitscore = calcBitScore(selectedHsp.score)
        eValue = calcEValue(len(seqHSP), dbSeqLength, bitscore)
        hspAlignment.dbStart += 1
        hspAlignment.seqStart += 2
        hspAlignment.seqEnd += 1
        temp = ""
        for j in range(len(seqHSP)):
            if j != 0:
                temp += seqHSP[j]
        seqHSP = temp

    posDB = hspAlignment.dbEnd
    if seqDB[posDB] != seqHSP[len(seqHSP)-1]:
        #On a un mismatch à la fin
        selectedHsp.score += 4
        bitscore = calcBitScore(selectedHsp.score)
        eValue = calcEValue(len(seqHSP), dbSeqLength, bitscore)
        hspAlignment.dbEnd -= 1
        hspAlignment.seqEnd -= 1

    print ("# Best HSP:")
    print ("Id:" + result.description + " Score brut:" + str(selectedHsp.score) + " Bitscore:" + str(
        bitscore) + " Evalue: " + str(eValue))
    dbStart = fixStartIndices(str(hspAlignment.dbStart))
    dbEnd = fixEndIndices(str(hspAlignment.dbEnd))
    inputStart = fixStartIndices(str(hspAlignment.seqStart))
    inputEnd = fixEndIndices(str(hspAlignment.seqEnd))
    print ("Alignement:")
    print (inputStart + seqHSP + inputEnd)
    print (dbStart + seqDB[hspAlignment.dbStart:(hspAlignment.dbEnd + 1)] + dbEnd)
    print ("---------------------------------------------")


def fixStartIndices(str):
    while len(str) < 3:
        str += " "
    return str


def fixEndIndices(str):
    while len(str) < 3:
        str = " " + str
    return str


def main():
    first = time.clock()
    parser = argparse.ArgumentParser(description='TP3_Bio')
    parser.add_argument('-i', help='Sequence ')
    parser.add_argument('-db', help='Base de donnees')
    parser.add_argument('-E', help='E Value')
    parser.add_argument('-ss', help='SS')
    parser.add_argument('-seed', '--seed', help='Graine 1 pour présent, 0 pour no care')
    args = parser.parse_args()
    if args.i:
        seqInput = args.i
        test = 0
        temp = ""
        for c in seqInput:
            if test != 0 and test != len(seqInput) - 1:
                temp += c
            test += 1
    else:
        return None
    seqInput = temp
    if args.db:
        seqSearchDB = fastaSequences(args.db)
    else:
        return None
    if args.E:
        seuil = float(args.E)
    else:
        seuil = 4
    if args.ss:
        seuilCutoff = float(args.ss)
    else:
        seuilCutoff = 0.001
    if args.seed:
        temp = ""
        test = 0
        for c in args.seed:
            if test != 0 and test != len(args.seed) - 1:
                temp += c
            test += 1
        seed = temp
    else:
        seed = "11111111111"

    selectedHspList = []
    kmerList = buildKmer(len(seed), seqInput)
    temp = ""
    for seqDB in seqSearchDB:
        if '>' in seqDB:
            temp = seqDB
            continue
        hspList, hspPos = findHSP(kmerList, seqDB, seed)
        hspExtendedList = extendGlouton(hspList, seqDB, seqInput, hspPos, seuil)
        hspMergedList = merge(hspExtendedList, seqInput, seqDB)
        dbSeqLength = getLengthSeqDB(seqSearchDB)
        selectedHsp = filterHSP(hspMergedList, len(seqInput), dbSeqLength, seuilCutoff)
        if selectedHsp is not None:
            selectedHspList.append(Result(selectedHsp, seqDB, temp))
    selectedHspList = quicksort(selectedHspList)

    for result in selectedHspList:
        selectedHsp = result.selectedHsp
        bitscore = selectedHsp.bitscore
        eValue = selectedHsp.eValue
        seqAlign = alignment(seqInput, result.seqDB)
        print (result.description + " Score: " + str(seqAlign.score) + " Ident: " + str(seqAlign.ident))
        printSmithWaterman(seqAlign)
        printAlignment(result, selectedHsp, bitscore, eValue, dbSeqLength)
    print ("Total : " + str(len(selectedHspList)))
    end = time.clock()
    print ("Time elapsed", end - first)



if __name__ == "__main__":
    main()
