import argparse
import os
import re
import sys
import unittest

match = 2
mismatch = -1
indel = -1

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

def main():
    seq2 = "GTCAGCC"
    seq1 = "CATAGTG"
    testAlignment = alignment(seq1, seq2)

if __name__ == "__main__":
    main()