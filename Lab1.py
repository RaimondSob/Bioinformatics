from typing import Sequence
from Bio import SeqIO
from Bio.Seq import Seq,reverse_complement, transcribe, back_transcribe, translate
from Bio.SeqRecord import SeqRecord



def findORF(seka,x):
     ORF = []
     for letter in range(x, len(seka), 3):
          endCodon = seka[letter:letter+3]
          if(endCodon == "TAA" or endCodon == "TAG" or endCodon == "TGA" ):
                    for letter1 in range(letter + 3, len(seka),3):
                         StopCodon = seka[letter1:letter1+3]
                         if(StopCodon == "TAA" or StopCodon == "TAG" or StopCodon == "TGA"): 
                              coding = seka[letter:letter1+3]
                              letter = letter1
                              if(len(coding)> 100):
                                   ORF.append(coding)
                              break
     return ORF

def findSeq(ORFS,x):
     CodingList = []
     for ORF in ORFS:
          for StartCodon in range(x, len(ORF), 3):
               cod = ORF[StartCodon:StartCodon+3]
               if(cod == "ATG"):
                    coding = ORF[StartCodon:len(ORF)]
                    #print(coding)
                    if(len(coding) > 100):
                         CodingList.append(coding)   
                    break
     return CodingList

def CountCodonFreq(seka):
     totalCodons = len(seka)/3
     Letters = ['A','T','C','G']
     CodonList = []
     for l1 in Letters:
          for l2 in Letters:
               for l3 in Letters:
                    codon = l1 + l2 + l3
                    CodonList.append(codon)

     CodonDict = dict.fromkeys(CodonList, 0.0)
     for letter in range(0, len(seka), 3):
          Codon = seka[letter: letter + 3]
          CodonDict[Codon] += 1
     for Codon in CodonDict:
          CodonDict[Codon] = CodonDict[Codon] / totalCodons * 100
     return CodonDict

def CountDicodonFreq(seka):
     if(len(seka) % 6 != 0):
          seka = seka[:-3]
     totalDicodons = len(seka)/6
     Letters = ['A','T','C','G']
     DicodonList = []
     for l1 in Letters:
          for l2 in Letters:
               for l3 in Letters:
                    for l4 in Letters:
                         for l5 in Letters:
                              for l6 in Letters:
                                   dicodon = l1 + l2 + l3 + l4 + l5 + l6
                                   DicodonList.append(dicodon) 
     DicodonDict = dict.fromkeys(DicodonList, 0.0)
     for letter in range(0, len(seka), 6):
          Dicodon = seka[letter: letter + 6]
          DicodonDict[Dicodon] += 1
     for Dicodon in DicodonDict:
          DicodonDict[Dicodon] = DicodonDict[Dicodon] / totalDicodons * 100
     return DicodonDict

def start_stop(seka,x):
     CodingList = []
     for letter in range(x, len(seka),3):
          startCodon = seka[letter:letter+3]
          if(startCodon == "ATG"):
               for letter1 in range(letter, len(seka),3):
                    endCodon = seka[letter1:letter1+3]
                    if(endCodon == "TAA" or endCodon == "TAG" or endCodon == "TGA"):
                         coding = seka[letter:letter1+3]
                         #print(coding)
                         letter = letter1 + 3 
                         if(len(coding) > 100):
                              CodingList.append(coding)
                         break
     return CodingList

def findCoding(Fasta):
     CodingList1 = findSeq(findORF(Fasta.seq,0), 0)
     CodingList2 = findSeq(findORF(Fasta.seq,1), 0)
     CodingList3 = findSeq(findORF(Fasta.seq,2), 0)
     Fasta2 = Seq(Fasta.seq.reverse_complement())
     CodingList4 = findSeq(findORF(Fasta2,0), 0)
     CodingList5 = findSeq(findORF(Fasta2,1), 0)
     CodingList6 = findSeq(findORF(Fasta2,2), 0)

     Coding1 = Seq('').join(CodingList1)
     Coding2 = Seq('').join(CodingList2)
     Coding3 = Seq('').join(CodingList3)
     Coding4 = Seq('').join(CodingList4)
     Coding5 = Seq('').join(CodingList5)
     Coding6 = Seq('').join(CodingList6)
     AllCoding = Coding1 + Coding2 + Coding3 + Coding4 + Coding5 + Coding6
     return AllCoding

Lactococcus_phage = SeqIO.read("C:\\Users\\Raimond\\Desktop\\viruses\\data\\bacterial1.fasta","fasta")
Escherichia_phage = SeqIO.read("C:\\Users\\Raimond\\Desktop\\viruses\\data\\bacterial2.fasta","fasta")
Streptococcus_phage = SeqIO.read("C:\\Users\\Raimond\\Desktop\\viruses\\data\\bacterial3.fasta","fasta")
Cellulophaga_phage = SeqIO.read("C:\\Users\\Raimond\\Desktop\\viruses\\data\\bacterial4.fasta","fasta")
coronavirus = SeqIO.read("C:\\Users\\Raimond\\Desktop\\viruses\\data\\mamalian1.fasta","fasta")
adenovirus = SeqIO.read("C:\\Users\\Raimond\\Desktop\\viruses\\data\\mamalian2.fasta","fasta")
Variola_virus = SeqIO.read("C:\\Users\\Raimond\\Desktop\\viruses\\data\\mamalian3.fasta","fasta")
herpesvirus = SeqIO.read("C:\\Users\\Raimond\\Desktop\\viruses\\data\\mamalian4.fasta","fasta")

Codon = CountCodonFreq(findCoding(Lactococcus_phage))

CodonFreq = []
CodonFreq.append(CountCodonFreq(findCoding(Lactococcus_phage)))
CodonFreq.append(CountCodonFreq(findCoding(Escherichia_phage)))
CodonFreq.append(CountCodonFreq(findCoding(Streptococcus_phage)))
CodonFreq.append(CountCodonFreq(findCoding(Cellulophaga_phage)))
CodonFreq.append(CountCodonFreq(findCoding(coronavirus)))
CodonFreq.append(CountCodonFreq(findCoding(adenovirus)))
CodonFreq.append(CountCodonFreq(findCoding(Variola_virus)))
CodonFreq.append(CountCodonFreq(findCoding(herpesvirus)))


Freq = []

CodonMatrix = [[0 for i in range(8)] for j in range (8)]
for i in range(0,8):
     for j in range(0,8):
          for key in Codon:
               Freq.append(abs(CodonFreq[i][key] - CodonFreq[j][key])/2)
          CodonMatrix[i][j] = sum(Freq)
          Freq = []
print(CodonMatrix)


Dicodon = CountDicodonFreq(findCoding(Lactococcus_phage))
DicodonFreq = []
DicodonFreq.append(CountDicodonFreq(findCoding(Lactococcus_phage)))
DicodonFreq.append(CountDicodonFreq(findCoding(Escherichia_phage)))
DicodonFreq.append(CountDicodonFreq(findCoding(Streptococcus_phage)))
DicodonFreq.append(CountDicodonFreq(findCoding(Cellulophaga_phage)))
DicodonFreq.append(CountDicodonFreq(findCoding(coronavirus)))
DicodonFreq.append(CountDicodonFreq(findCoding(adenovirus)))
DicodonFreq.append(CountDicodonFreq(findCoding(Variola_virus)))
DicodonFreq.append(CountDicodonFreq(findCoding(herpesvirus)))

Freq = []

DicodonMatrix = [[0 for i in range(8)] for j in range (8)]
for i in range(0,8):
     for j in range(0,8):
          for key in Dicodon:
               Freq.append(abs(DicodonFreq[i][key] - DicodonFreq[j][key])/2)
          DicodonMatrix[i][j] = sum(Freq)
          Freq = []
print(DicodonMatrix)
