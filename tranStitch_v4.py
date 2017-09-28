stats = {}
import sys
def extractScaffolds(assembly,reassemblyFile):
    seqDict = buildGenomeDict(assembly)
    infile = open(reassemblyFile,'r')
    scaffoldList = []
    for line in infile:
        lineSplit = line.split('\t')
        scaff1 = '>' + lineSplit[0]
        scaff2 = '>' + lineSplit[1]
        if scaff1 not in scaffoldList:
            scaffoldList.append(scaff1)
        if scaff2 not in scaffoldList:
            scaffoldList.append(scaff2)
    infile.close()
    outfile = open('scaffoldsWithDuplications_v3_' + assembly[0:-6] + '.fasta','w')
    stats['# Reassembly Recs'] = len(scaffoldList)
    for seq in scaffoldList:
        outfile.write(seq + '\n')
        outfile.write(seqDict[seq] + '\n')
    outfile.close()
    return stats
    
def reverseComplement(seq):
    seq_revc = ''
    for nuc in seq:
        if nuc == 'A':
            seq_revc = 'T' + seq_revc
        elif nuc == 'T':
            seq_revc = 'A' + seq_revc
        elif nuc == 'C':
            seq_revc = 'G' + seq_revc
        elif nuc == 'G':
            seq_revc = 'C' + seq_revc
        else:
            seq_revc = nuc + seq_revc
    return seq_revc
    
def buildSeqDict(listOfAssemblyFiles):
    seqDict = {}
    seqName = ''
    currSeq = ''
    for assembly in listOfAssemblyFiles:
        infile = open(assembly,'r')
        fileNameSplit = assembly.split('_')
        assemblyTag = fileNameSplit[0]
        for line in infile:
            if line[0] == '>':
                if seqName != '':
                    seqDict[seqName] = currSeq
                seqName = line
                while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                    seqName = seqName[0:-1]
                seqName  = assemblyTag + '_' + seqName[1:]
                currSeq = ''
            else:
                currSeq += line.upper()
                while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                    currSeq = currSeq[0:-1]
        seqDict[seqName] = currSeq
        infile.close()
    return seqDict

def buildGenomeDict(assembly):
    seqDict = {}
    seqName = ''
    currSeq = ''
    infile = open(assembly,'r')
    for line in infile:
        if line[0] == '>':
            if seqName != '':
                seqDict[seqName] = currSeq
            seqName = line
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
            currSeq = ''
        else:
            currSeq += line.upper()
            while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                currSeq = currSeq[0:-1]
    seqDict[seqName] = currSeq
    infile.close()
    return seqDict
    
def orient(lineSplit, scaffoldDict, transcriptDict):              
    logfile = open('logfile.txt','a')
    scaffoldTag1 = '>' + lineSplit[1]
    scaffoldTag2 = '>' + lineSplit[6]
    scaffold1 = scaffoldDict[scaffoldTag1]
    scaffold2 = scaffoldDict[scaffoldTag2]
    orient1 = int(lineSplit[3])
    orient2 = int(lineSplit[8])
    transcriptTag = lineSplit[0]
    counter = 0
    evidence = []
    for item in lineSplit[0:-1]:
        if counter%11 == 0:
            evidence.append(item)
        counter += 1
    currTranscript = transcriptDict[transcriptTag]
    alignment1 = lineSplit[2]
    alignment2 = lineSplit[7]
    if orient1 == 16 or orient1 == 272:
        alignment1 = reverseComplement(alignment1)
        alignment1 = alignment1[-12:]
    else:
        alignment1 = alignment1[0:12]
    if orient2 == 16 or orient2 == 272:
        alignment2 = reverseComplement(alignment2)
        alignment2 = alignment2[-12:]
    else:
        alignment2 = alignment2[0:12]
    numAlignments1 = 0
    numAlignments2 = 0
    for i in range(len(currTranscript)-11):
        currWord = currTranscript[i:i+12]
        if alignment1 == currWord:
            transPos1 = i
            numAlignments1 += 1
        if alignment2 == currWord:
            transPos2 = i
            numAlignments2 += 1
    if orient1 == 16 or orient1 == 272:
        transPos1 = len(currTranscript) - transPos1 - 12
    if orient2 == 16 or orient2 == 272:
        transPos2 = len(currTranscript) - transPos2 - 12
    if transPos1 < transPos2:
        if orient1 == orient2:
            newScaffold = scaffold1 + 'N'*50 + scaffold2
            newScaffoldTag = '>' + scaffoldTag1[1:] + '_' + scaffoldTag2[1:] + '_stitched'
        elif orient1 == 0 or orient1 == 256:
            if orient2 == 0 or orient2 == 256:
                newScaffold = scaffold1 + 'N'*50 + scaffold2
                newScaffoldTag = '>' + scaffoldTag1[1:] + '_' + scaffoldTag2[1:] + '_stitched'
            else:
                scaffold2_revc = reverseComplement(scaffold2)
                newScaffold = scaffold1 + 'N'*50 + scaffold2_revc
                newScaffoldTag = '>' + scaffoldTag1[1:] + '_' + scaffoldTag2[1:] + '_rev_comp_stitched'
        else:
            if orient2 == 16 or orient2 == 272:
                newScaffold = scaffold1 + 'N'*50 + scaffold2
                newScaffoldTag = '>' + scaffoldTag1[1:] + '_' + scaffoldTag2[1:] + '_stitched'
            else:
                scaffold1_revc = reverseComplement(scaffold1)
                newScaffold = scaffold1_revc + 'N'*50 + scaffold2
                newScaffoldTag = '>' + scaffoldTag1[1:] + '_rev_comp_' + scaffoldTag2[1:] + '_stitched'
        logfile.write('Stitched scaffolds ' + scaffoldTag2[1:] + ' and ' + scaffoldTag1[1:] + ' based on the following transcripts:' + '\n')
        logfile.write('\t')
        for item in evidence:
            logfile.write(item + '\t')
        logfile.write('\n')
    elif transPos1 > transPos2:
        if orient1 == orient2:
            newScaffold = scaffold2 + 'N'*50 + scaffold1
            newScaffoldTag = '>' + scaffoldTag2[1:] + '_' + scaffoldTag1[1:] + '_stitched'
        elif orient1 == 0 or orient1 == 256:
            if orient2 == 0 or orient2 == 256:
                newScaffold = scaffold2 + 'N'*50 + scaffold1
                newScaffoldTag = '>' + scaffoldTag2[1:] + '_' + scaffoldTag1[1:] + '_stitched'
            else:
                scaffold2_revc = reverseComplement(scaffold2)
                newScaffold = scaffold2_revc + 'N'*50 + scaffold1
                newScaffoldTag = '>' + scaffoldTag2[1:] + '_rev_comp_' + scaffoldTag1[1:] + '_stitched'
        else:
            if orient2 == 16 or orient2 == 272:
                newScaffold = scaffold2 + 'N'*50 + scaffold1
                newScaffoldTag = '>' + scaffoldTag2[1:] + '_' + scaffoldTag1[1:] + '_stitched'
            else:
                scaffold1_revc = reverseComplement(scaffold1)
                newScaffold = scaffold1_revc + 'N'*50 + scaffold2
                newScaffoldTag ='>' + scaffoldTag2[1:] + '_' + scaffoldTag1[1:] + '_rev_comp_stitched'
        logfile.write('Stitched scaffolds ' + scaffoldTag2[1:] + ' and ' + scaffoldTag1[1:] + ' based on the following transcripts:\n')
        logfile.write('\t')
        for item in evidence:
            logfile.write(item + '\t')
        logfile.write('\n')
        logfile.close()
        return newScaffoldTag,newScaffold
    else:
        logfile.write('Could not determine stitching orientation for scaffolds ' + scaffoldTag1 + ' and ' + scaffoldTag2 + '\n')
        logfile.close()
        newScaffoldTag,newScaffold = 'No stitch','No stitch'
    return newScaffoldTag,newScaffold

def mergeStitches(stitchFile, genome, listOfTranscriptomes):
    stitchRecs = open(stitchFile,'r')
    scaffoldDict = buildGenomeDict(genome)
    transcriptDict = buildSeqDict(listOfTranscriptomes)
    logfile = open('logfile.txt','a')
    logfile.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n' + genome[0:-6] + '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    scaffoldList = []
    stitchDict = {}
    for line in stitchRecs:
        lineSplit = line.split('\t')
        scaffoldTag1 = '>' + lineSplit[1]
        scaffoldTag2 = '>' + lineSplit[6]
        if scaffoldTag1 not in scaffoldList:
            scaffoldList.append(scaffoldTag1)
        if scaffoldTag2 not in scaffoldList:
            scaffoldList.append(scaffoldTag2)
        numTranscripts = len(lineSplit)/11
        stitchDict[(scaffoldTag1,scaffoldTag2)] = lineSplit + [numTranscripts]
    stitchRecs.close()
    stats['# Possible Stitches'] = len(stitchDict)
    numStitches = 0
    for a in scaffoldList:
        currPairs = []
        for b in scaffoldList:
            if (a,b) in stitchDict:
                stitchInfo = stitchDict[(a,b)]
                currPairs.append([a,b,stitchInfo[-1]])
            elif (b,a) in stitchDict:
                stitchInfo = stitchDict[(b,a)]
                currPairs.append([b,a,stitchInfo[-1]])
        if len(currPairs) == 1:
            pair = currPairs[0]
            stitchInfo = stitchDict[(pair[0],pair[1])]
            currStitchName, currStitch = orient(stitchInfo, scaffoldDict, transcriptDict)
            if currStitchName != 'No stitch':
                del scaffoldDict[pair[0]]
                del scaffoldDict[pair[1]]
                scaffoldList.remove(pair[0])
                scaffoldList.remove(pair[1])
                scaffoldDict[currStitchName] = currStitch
                numStitches += 1
        else:
            maxEvidence = 0
            sortedIndices = []
            index = 0
            currPairIndex = -1
            for pair in currPairs:
                currEvidence = pair[2]
                if currEvidence >= maxEvidence:
                    maxEvidence = currEvidence
                    sortedIndices.append(index)
                else:
                    i = 0
                    for item in sortedIndices:
                        pair = currPairs[item]
                        numTranscripts = pair[2]
                        if currEvidence <= numTranscripts:
                            sortedIndices = sortedIndices[0:i] + [index] + sortedIndices[i:]
                        elif i < len(sortedIndices) - 1:
                            i += 1
                        else:
                            sortedIndices.append(index)
                index += 1
            sortedPairs = []
            for num in sortedIndices:
                pair = currPairs[num]
                sortedPairs.append(pair)
            newPairs = []
            for pair in currPairs:
                if pair[2] == maxEvidence:
                    newPairs.append(pair)
            if len(newPairs) > 1:
                stitchLengths = []
                maxLength = 0
                for pair in newPairs:
                    stitchLength = len(scaffoldDict[pair[0]]) + len(scaffoldDict[pair[1]])
                    if stitchLength > maxLength:
                        stitchLengths.append(pair[0:2] + [stitchLength])
                        maxLength = stitchLength
                stitches = []
                for stitch in stitchLengths:
                    if stitch[2] == maxLength:
                        stitches.append((stitch[0],stitch[1]))
                stitchInfo = stitchDict[stitches[0]]
                currStitchName, currStitch = orient(stitchInfo, scaffoldDict, transcriptDict)
                if currStitchName != 'No stitch':
                    del scaffoldDict[pair[0]]
                    del scaffoldDict[pair[1]]
                    scaffoldList.remove(pair[0])
                    scaffoldList.remove(pair[1])
                    scaffoldDict[currStitchName] = currStitch
                    numStitches += 1
                else:
                    while currStitchName == 'No stitch' and -currPairIndex <  len(currPairs):
                        currPairIndex -= 1
                        nextStitch = sortedPairs[currPairIndex]
                        stitchInfo = stitchDict[(nextStitch[0],nextStitch[1])]
                        currStitchName, currStitch = orient(stitchInfo, scaffoldDict, transcriptDict)
                    if currStitchName != 'No stitch':
                        del scaffoldDict[pair[0]]
                        del scaffoldDict[pair[1]]
                        scaffoldList.remove(pair[0])
                        scaffoldList.remove(pair[1])
                        scaffoldDict[currStitchName] = currStitch
                        numStitches += 1    
            elif len(newPairs) == 1:
                pair = newPairs[0]
                stitchInfo = stitchDict[(pair[0],pair[1])]
                currStitchName, currStitch = orient(stitchInfo, scaffoldDict, transcriptDict)
                if currStitchName != 'No stitch':
                    del scaffoldDict[pair[0]]
                    del scaffoldDict[pair[1]]
                    scaffoldList.remove(pair[0])
                    scaffoldList.remove(pair[1])
                    scaffoldDict[currStitchName] = currStitch
                    numStitches += 1
                else:
                    while currStitchName == 'No stitch' and (-currPairIndex) <  len(currPairs):
                        currPairIndex -= 1
                        nextStitch = sortedPairs[currPairIndex]
                        stitchInfo = stitchDict[(nextStitch[0],nextStitch[1])]
                        currStitchName, currStitch = orient(stitchInfo, scaffoldDict, transcriptDict)
                    if currStitchName != 'No stitch':
                        del scaffoldDict[pair[0]]
                        del scaffoldDict[pair[1]]
                        scaffoldList.remove(pair[0])
                        scaffoldList.remove(pair[1])
                        scaffoldDict[currStitchName] = currStitch 
                        numStitches += 1   
    outfile = open('stitched_scaffolds_v3_' + genome[0:-6] + '.fasta','w')
    lengthList = []
    lengthDict = {}
    for seq in scaffoldDict:
        currLength = len(scaffoldDict[seq])
        if currLength not in lengthList:
            lengths = [seq]
            lengthDict[currLength] = lengths
            lengthList.append(currLength)
        else:
            lengths = lengthDict[currLength]
            lengths.append(seq)
            lengthDict[currLength] = lengths
    lengthList = sorted(lengthList, reverse = True)
    for length in lengthList:
        seqs = lengthDict[length]
        for seq in seqs:
            outfile.write(seq + '\n')
            outfile.write(scaffoldDict[seq] + '\n')
    outfile.close()
    stats['# Stitches Made'] = numStitches
    
def compareAlignments(seq1,seq2):
    same = True
    if len(seq1) <= len(seq2):
        for i in range(len(seq1)): 
            if seq1[i] != seq2[i]:
                same = False
        if same == False:
            seq2_revc = ''
            for nuc in seq2:
                if nuc == 'A':
                    seq2_revc = 'T' + seq2_revc
                elif nuc == 'T':
                    seq2_revc = 'A' + seq2_revc
                elif nuc == 'C':
                    seq2_revc = 'G' + seq2_revc
                elif nuc == 'G':
                    seq2_revc = 'C' + seq2_revc
                else:
                    seq2_revc = nuc + seq2_revc
            same = True
            for i in range(len(seq1)):
                if seq1[i] != seq2_revc[i]:
                    same = False
    else:
        for i in range(len(seq2)): 
            if seq1[i] != seq2[i]:
                same = False
        if same == False:
            seq2_revc = ''
            for nuc in seq2:
                if nuc == 'A':
                    seq2_revc = 'T' + seq2_revc
                elif nuc == 'T':
                    seq2_revc = 'A' + seq2_revc
                elif nuc == 'C':
                    seq2_revc = 'G' + seq2_revc
                elif nuc == 'G':
                    seq2_revc = 'C' + seq2_revc
                else:
                    seq2_revc = nuc + seq2_revc
            same = True
            for i in range(len(seq2)):
                if seq1[i] != seq2_revc[i]:
                    same = False
    return same
            

def stitch(genome, listOfSamFiles):
    outfile = open('stitch_v3_Recommendations_' + genome[0:-6] + '.txt','w')
    outfile2 = open('duplicatedRegions_v3_' + genome[0:-6] + '.txt','w')
    scaffoldDict = {}
    scaffoldList = []
    for infile in listOfSamFiles:
        currFile = open(infile,'r')
        nameSplit = infile[0:-4].split('_')
        for line in currFile:
            if line[0] != '@':
                lineSplit = line.split('\t')
                transcriptTag = nameSplit[-1] + '_' + lineSplit[0]#Critical for proper stitching!
                transcriptSplit = transcriptTag.split('_')
                if len(transcriptSplit) > 4:
                    transcriptTag = transcriptSplit[0] + '_'
                    for i in range(1,len(transcriptSplit)-1):
                        transcriptTag += transcriptSplit[i] + '_'
                    transcriptTag = transcriptTag[0:-1]
                orientation = lineSplit[1]
                scaffold = lineSplit[2]
                startPos = lineSplit[3]
            	numMapLocationTag = int(lineSplit[4])
                cigar = lineSplit[5]
                alignment = lineSplit[9]    
                if numMapLocationTag == 50 or numMapLocationTag == 3: #or numMapLocationTag == 3: #remove second condition for real code
                    transcriptInfo = [alignment,orientation,startPos,cigar]
                    scaffIn = False
                    for scaff in scaffoldList:
                        if scaff == scaffold:
                            scaffIn = True
                    if scaffIn == True:
                        currDict = scaffoldDict[scaffold]
                        currDict[transcriptTag] = transcriptInfo
                        scaffoldDict[scaffold] = currDict
                    else:
                        currDict = {transcriptTag:transcriptInfo}
                        scaffoldDict[scaffold] = currDict
                        scaffoldList.append(scaffold)
        currFile.close()
    stitch = False
    counter = 1
    stitchDict = {}
    stitchScaffoldList = []
    numScaffoldsSharingTranscript = 0
    for scaffold1 in scaffoldList[0:-1]:
        scaff1Dict = scaffoldDict[scaffold1]
        scaff1TranscriptList = []
        for key in scaff1Dict:
            scaff1TranscriptList.append(key)
        for scaffold2 in scaffoldList[counter:]:
            scaff2Dict = scaffoldDict[scaffold2]
            scaff2TranscriptList = []
            for key in scaff2Dict:
                scaff2TranscriptList.append(key)
            stitchList = []
            for trans1 in scaff1TranscriptList:
                transcript1Info = scaff1Dict[trans1]
                for trans2 in scaff2TranscriptList:
                    stitch = False
                    transcript1Info = scaff1Dict[trans1]
                    align1 = transcript1Info[0]
                    transcript2Info = scaff2Dict[trans2]
	            align2 = transcript2Info[0]
                    if trans1 == trans2:
                        test = compareAlignments(align1,align2)
                        numScaffoldsSharingTranscript += 1
                        if test == False:
                            stitch = True
                    if stitch == True:
                        stitchList.append(trans1)
                        stitchDict[(scaffold1,scaffold2)] = stitchList
                        stitchScaffoldList.append((scaffold1,scaffold2))
                    elif trans1 == trans2:
                        startPos1 = int(transcript1Info[2])
                        cigar1 = transcript1Info[3]
                        startPos2 = int(transcript2Info[2])
                        cigar2 = transcript2Info[3]
                        orient1 = int(transcript1Info[1])
                        orient2 = int(transcript2Info[1])
                        totalLength = 0
                        currString = ''
                        for char in cigar1:
                            if char == 'M':
                                totalLength += int(currString)
                                currString = ''
                            elif char == 'N':
                                totalLength += int(currString)
                                currString = ''
                            elif char == 'I':
                                totalLength += int(currString)
                                currString = ''
                            elif char == 'D':
                                currString = ''
                            elif char == 'S':
                                totalLength += int(currString)
                                currString = ''
                            elif char == 'H':
                                currString = ''
                            elif char == 'P':
                                currString = ''
                            elif char == 'X':
                                totalLength += int(currString)
                                currString = ''
                            elif char == '=':
                                totalLength += int(currString)
                                currString = ''
                            else:
                                currString += char
                        if orient1 == 16 or orient1 == 272:
                            endPos1 = startPos1 - totalLength
                        elif orient1 == 0 or orient1 == 256:
                            endPos1 = startPos1 + totalLength
                        else:
                            endPos1 = False
                        totalLength = 0
                        currString = ''
                        for char in cigar2:
                            if char == 'M':
                                totalLength += int(currString)
                                currString = ''
                            elif char == 'N':
                                totalLength += int(currString)
                                currString = ''
                            elif char == 'I':
                                totalLength += int(currString)
                                currString = ''
                            elif char == 'D':
                                currString = ''
                            elif char == 'S':
                                totalLength += int(currString)
                                currString = ''
                            elif char == 'H':
                                currString = ''
                            elif char == 'P':
                                currString = ''
                            elif char == 'X':
                                totalLength += int(currString)
                                currString = ''
                            elif char == '=':
                                totalLength += int(currString)
                                currString = ''
                            else:
                                currString += char
                        if orient2 == 16 or orient2 == 272:
                            endPos2 = startPos2 - totalLength
                        elif orient2 == 0 or orient2 == 256:
                            endPos2 = startPos2 + totalLength
                        else:
                            endPos2 = False
                        if endPos1 != False and endPos2 != False:
                            outfile2.write(trans1 + '\t' + scaffold1 + '\t' + str(orient1) + '\t' + str(startPos1) + '\t' + str(endPos1) + '\t' + scaffold2 + '\t' + str(orient2)  + '\t' + str(startPos2) + '\t' + str(endPos2) + '\n')
        counter += 1
    for stitch in stitchScaffoldList:
        scaff1 = stitch[0]
        scaff2 = stitch[1]
        scaff1Dict = scaffoldDict[scaff1]
        scaff2Dict = scaffoldDict[scaff2]
        evidence = stitchDict[stitch]
        for transcript in evidence:
            transcriptInfo1 = scaff1Dict[transcript]
            transcriptInfo2 = scaff2Dict[transcript]
            outfile.write(transcript + '\t')
            outfile.write(scaff1 + '\t')
            for item in transcriptInfo1:
                outfile.write(item + '\t')
            outfile.write(scaff2 + '\t')
            for item in transcriptInfo2:
                outfile.write(item + '\t')
        outfile.write('\n')
    stats['# Scaffolds Mapped To'] = len(scaffoldList)
    stats['# Scaffolds Sharing a Transcript'] = numScaffoldsSharingTranscript
    outfile.close()
    outfile2.close()
    return stats
                        
                    
                        
                              
    
def transcriptMultiMap(samFile):
    infile = open(samFile,'r')
    transcriptDict = {}
    transcriptList = []
    scaffoldDict = {}    
    for line in infile:
        if line[0] != '@':
            lineSplit = line.split('\t')
            transcriptTag = lineSplit[0]
            transcriptSplit = transcriptTag.split('_')
            if len(transcriptSplit) > 3:
                transcriptTag = transcriptSplit[0]
                for i in range(1,len(transcriptSplit)-1):
                    transcriptTag += transcriptSplit[i] + '_'
                transcriptTag = transcriptTag[0:-1]
            orientation = lineSplit[1]
            scaffold = lineSplit[2]
            startPos = lineSplit[3]
            numMapLocationTag = int(lineSplit[4])
            cigar = lineSplit[5]
            if numMapLocationTag == 50 or numMapLocationTag == 3:
                if transcriptTag in transcriptList:
                    currDict = transcriptDict[transcriptTag]
                    scaffIn = False
                    i = 0
                    while scaffIn == False and i < len(currDict):
                        for scaff in currDict:
                            if scaff == scaffold:
                                scaffIn == True
                            i += 1
                    if scaffIn == False:
                        scaffInfo = [[orientation, startPos, cigar]]
                        currDict[scaffold] = scaffInfo
                    else:
                        scaffInfo = [[orientation, startPos, cigar]]
                        currDict[scaffold] += scaffInfo
                    
                    transcriptDict[transcriptTag] = currDict
                else:
                    transcriptList.append(transcriptTag)
                    currDict = {}
                    scaffInfo = [[orientation, startPos, cigar]]
                    currDict[scaffold] = scaffInfo
                    transcriptDict[transcriptTag] = currDict
    outfile = open('transcript_supports_' + samFile[0:-4] + '.txt','w')
    for transcript in transcriptList:
        if len(transcriptDict[transcript]) > 1:
            scaffDict = transcriptDict[transcript]
            outfile.write(transcript + '\t')
            for scaffold in scaffDict:
                outfile.write(scaffold + '\t')
                for transcriptRegion in scaffDict[scaffold]:
                    for item in transcriptRegion:
                        outfile.write(str(item) + '\t')
            outfile.write('\n')
    infile.close()
    outfile.close()


def scaffolds2Stitch(Infile):
    infile = open(Infile,'r')
    outfile = open('scaffolds2stitch_' + Infile[20:-4] + '.txt','w')
    scaffoldDict = {}
    for line in infile:
        scaffoldPairs = []
        if line[0] != 'T':
            scaffoldList = []
            lineSplit = line.split('\t')
            i = 0
            scaffoldList.append(lineSplit[1])
            for item in lineSplit[1:-1]:
                if i%4 == 0:
                    scaffoldList.append(item)
                i += 1    
            for i in range(len(scaffoldList)):
                for j in range(len(scaffoldList)):
                    if scaffoldList[i] != scaffoldList[j]:
                        if (scaffoldList[i],scaffoldList[j]) not in scaffoldPairs and (scaffoldList[j],scaffoldList[i]) not in scaffoldPairs:
                            scaffoldPairs.append((scaffoldList[i],scaffoldList[j]))
        for pair in scaffoldPairs:
            if pair not in scaffoldDict:
                scaffoldDict[pair] = [lineSplit[0]]
            else:
                scaffoldDict[pair] += [lineSplit[0]]
    outfile.write('Scaffold 1' + '\t' + 'Scaffold 2' + '\t' + 'Transcripts Supporting' + '\n')
    for pair in scaffoldDict:
        print pair
        outfile.write(pair[0] + '\t' + pair[1])
        transcriptList = scaffoldDict[pair]
        for transcript in transcriptList:
            outfile.write('\t' + transcript)
        outfile.write('\n')
    infile.close()
    outfile.close()

import sys
def matchMisMatch(i,j):
    if i == j:
        return (2)
    else:
        return (-1)

def align(Infile):
    infile = open(Infile,'r')
    outfile = open('alignment.fasta','w')
    seqDict = {}
    currentSeq = ''
    currentSeqName = ''
    counter = 0
    for line in infile:
        counter += 1
        if line[0] == '>' and currentSeqName == '':
            currentSeqName = line[0:-1]
        elif line[0] == '>':
            seqDict[currentSeqName] = currentSeq
            currentSeqName = line
            while currentSeqName[-1] == '\n' or currentSeqName[-1] == '\t' or currentSeqName[-1] == '\r':
                currentSeqName = currentSeqName[0:-1]           
            currentSeq = ''
        else:
            currentSeq += line
            while currentSeq[-1] == '\n' or currentSeq[-1] == '\t' or currentSeq[-1] == '\r':
                currentSeq = currentSeq[0:-1]
            currentSeq = currentSeq.upper()
    seqDict[currentSeqName] = currentSeq
    seq1 = ''
    for i in seqDict:
        if seq1 == '':
            seq1 = seqDict[i]
        else:
            seq2 = seqDict[i]
    a = {}
    m = len(seq1)
    n = len(seq2)
    gapScore = (-2)
    gapExtend = 1
    for i in range(m+1):
        a[i,0] = i*gapScore
    for j in range(n+1):
        a[0,j] = j*gapScore
    for i in range(1,m + 1):
        for j in range(1,n + 1):
            a[i,j] = max((a[i-1,j] + gapScore),(a[i-1,j-1] + matchMisMatch(seq1[i-1],seq2[j-1])),(a[i,j-1] + gapScore))
    x = m
    y = n                 
    newSeq1 = ''
    newSeq2 = ''
    while x > 0 and y > 0:
         threeSquares = [a[x,y-1],a[x-1,y-1],a[x-1,y]]
         if threeSquares[0] > threeSquares[1] and threeSquares[0] > threeSquares[2]:
             newSeq1 = '-' + newSeq1
             newSeq2 = seq2[y-1] + newSeq2
             y = y - 1
         elif threeSquares[1] >= threeSquares[2]:
             newSeq1 = seq1[x-1] + newSeq1
             newSeq2 = seq2[y-1] + newSeq2
             x = x - 1
             y = y - 1
         else:
             newSeq1 = seq1[x-1] + newSeq1
             newSeq2 = '-' + newSeq2
             x = x - 1
    outfile.write('>Seq1\n')
    outfile.write(newSeq1 + '\n')
    outfile.write('>Seq2\n')
    outfile.write(newSeq2 + '\n')
    infile.close()
    outfile.close()
                    
def runCommand():
    if sys.argv[1] == 'stitch':
        sams = sys.argv[3]
        if sams[-4:] == 'fofn':
            samList = []
            infile = open(sams,'r')
            for line in infile:
                samList.append(line[0:-1])
            infile.close()
            stitch(sys.argv[2],samList)
        elif sams[-3:] == 'sam':
            samList = [sams]
            stitch(sys.argv[2],samList)
        else:
            print 'Bad formatting'
            return 'Bad formatting'
    elif sys.argv[1] == 'merge':
        reads = sys.argv[4]
        if reads[-4:] == 'fofn':
            readList = []
            infile = open(reads,'r')
            for line in infile:
                readList.append(line[0:-1])
            infile.close()
            mergeStitches(sys.argv[2],sys.argv[3],readList)
        elif reads[-5:] == 'fasta':
            mergeStitches(sys.argv[2],sys.argv[3],sys.argv[4])
        else:
            print 'Bad formatting'
            return 'Bad formatting'
    else:
        print 'No command specified'
        return 'No command specified'

runCommand()                
    
