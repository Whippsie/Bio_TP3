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
    for hsp in hspList:
        rightString = ""
        leftString = ""
        # TODO: Probleme si plus qu'un HSP dans la sequence, nous donnera toujours le meme
        posDB = hspList.index(hsp)
        posDB = hspPos[posDB]
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


# CODE TP1
def alignLocal(seq1, seqDB):
    matrix, scoreMax, scoreMaxPos = makeMatrix(seq1, seqDB)
    print(matrix)
    print(scoreMax)
    path, end = sequencePath(matrix, scoreMaxPos, seq1, seqDB)
    aligned, alignValue = alignSequences(end, path, seq1, seqDB, scoreMaxPos, matrix.shape)
    print(aligned)
    print("Alignement value:", alignValue)


def makeMatrix(seq1, seqDB):
    matrix = [[0 for col in range(len(seq1)+1)] for row in range(len(seqDB+1))]
    shape = matrix.shape
    indel = -1
    scoreMax = 0
    scoreMaxPos = (0, 0)
    for i in range(1, shape[0]):
        for j in range(1, shape[1]):
            score = getScore(seq1, seqDB, i, j, matrix)
            if checkMax(score, scoreMax):
                scoreMax = score
                scoreMaxPos = (i, j)
    return matrix, scoreMax, scoreMaxPos


def alignSequences(start, path, seq1, seq2, end, size):
    path.reverse()
    temp = seq2
    seq2 = seq1
    seq1 = temp

    """Indels de départ"""
    seq1align, seq2align, x, y, z = genIndelStart(start, seq1, seq2)

    """Séquences de chemin"""
    seq1align, seq2align, y, z, i = genSeqPath(seq1align, seq2align, path, y, z, seq1, seq2)

    """Complétion de la séquence"""
    seq1align, seq2align = genIndelEnd(seq1align, seq2align, end, size, seq1, seq2)

    return [seq1align, seq2align], i


def getScore(seq1, seqDB, i, j, matrix):
    if seq1[i - 1] == seqDB[j - 1]:
        mymatch = 1
    else:
        mymatch = -1
    dscore = matrix[i - 1][j - 1] + mymatch
    uscore = matrix[i - 1][j] + indel
    lscore = matrix[i][j - 1] + indel
    return max(0, dscore, uscore, lscore)


def checkMax(score, scoreMax):
    return score > scoreMax


def sequencePath(matrice, pos, seq1, seq2):
    x = pos[0]
    y = pos[1]
    path = []
    current = matrice[x][y]

    while ((x > 0) and (y > 0)):
        if ((current - mismatch == matrice[x - 1][y - 1]) and (seq1[x - 1] != seq2[y - 1])) or (
            (current - match == matrice[x - 1][y - 1]) and (seq1[x - 1] == seq2[y - 1])):
            x -= 1
            y -= 1
            path.append([1, 1])
        elif (current - indel == matrice[x - 1][y]):
            x -= 1
            current = matrice[x][y]
            path.append([1, 0])
        elif (current - indel == matrice[x][y - 1]):
            y -= 1
            current = matrice[x][y]
            path.append([0, 1])
        else:
            break
        current = matrice[x][y]
    #return path, numpy.array([x, y])
    return path, None

def genIndelStart(start, seq1, seq2):
    y = 0
    z = 0
    seq1align = ""
    seq2align = ""
    if start[0] > 0:
        x = start[0]
        while x > 0:
            seq1align += "-"
            x -= 1
            seq2align += seq2[y]
            y += 1
    if start[1] > 0:
        x = start[1]
        while x > 0:
            seq2align += "-"
            x -= 1
            seq1align += seq1[z]
            z += 1
    else:
        pass
    return seq1align, seq2align, x, y, z


def genSeqPath(s1, s2, p, y, z, sq1, sq2):
    i = 0
    while i < len(p):
        if p[i][0] == 1:
            s2 += sq2[y]
            y += 1
        else:
            s2 += "-"

        if p[i][1] == 1:
            s1 += sq1[z]
            z += 1
        else:
            s1 += "-"
        i += 1

    return s1, s2, y, z, i


def genIndelEnd(s1, s2, end, size, sq1, sq2):
    size_ligne = size[0] - 1
    size_col = size[1] - 1

    if end[0] < size_ligne:
        x = end[0]
        while x < size_ligne:
            s1 += "-"
            s2 += sq2[x]
            x += 1
    elif end[1] < size_col:
        x = end[1]
        while x < size_col:
            s2 += "-"
            s1 += sq1[x]
            x += 1

    return s1, s2


def reverseSeq(seq):
    seqrev = ""
    for s in seq:
        if s == "A":
            seqrev += "T"
        elif s == "C":
            seqrev += "G"
        elif s == "T":
            seqrev += "A"
        elif s == "G":
            seqrev += "C"

    l = list(seqrev)
    l.reverse()
    seqrev = ''.join(l)

    return seqrev


# CODE TP1

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
