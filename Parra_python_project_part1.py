#!/usr/bin/env python
# coding: utf-8

# # A. Make a class called seq
#
# ## This class should accept the following attributes:
# ### name, organism, sequence, and type. When making an isntance of seq, the input for name will be a string (e.g. 'RAS_G12D') the input for organism will be a string (e.x. 'human'), the input for sequence will be a DNA, RNA, or protein sequence (e.g. 'ATCGAAATC') and the type will be either 'DNA', 'RNA', or 'Protein'
#
# ## This class should have three methods:
# ### 1. info -this should print the name, type, organism, and sequence of the instance
# ### 2. length -this should count the length of the sequence string
# ### 3. fasta_out -this should write the name, organism, type, and sequence as a fasta file.

# In[137]:


# Helpful variables

standard_code = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "UGA": "*",
    "UGU": "C",
    "UGC": "C",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

aa_mol_weights = {
    "A": 89.09,
    "C": 121.15,
    "D": 133.1,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "K": 146.19,
    "L": 131.17,
    "M": 149.21,
    "N": 132.12,
    "P": 115.13,
    "Q": 146.15,
    "R": 174.2,
    "S": 105.09,
    "T": 119.12,
    "V": 117.15,
    "W": 204.23,
    "X": 0,
    "Y": 181.19,
}


# In[138]:


# call the class seq
class seq:

    # Call instance attributes
    def __init__(self, name, organism, sequence, type):
        self.name = name
        self.organism = organism
        self.sequence = sequence
        self.type = type

    # define the info function
    def info(self):
        print(self.name)
        print(self.type)
        print(self.organism)
        print(self.sequence)

    # define the length function
    def length(self):
        return len(self.sequence)

    # define the fasta_out function.
    # write to a file using the sequence name as part of the file name
    # I have written most of the code for the fasta_out function,
    # but you need to add something. What is missing?
    # this url may be helpful https://www.w3schools.com/python/python_file_write.asp
    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "\n"
            + self.sequence
            + "\n"
        )
        f.close()


# # C. Using super(), make a new child class of seq called protein
# ## This should have a new attribute called size (For instances of protein, size values will be in kDa, like '52')
# ## Overwrite the parent class seq function fasta_out to include the protein size in the first line of the fasta file
#

# In[140]:


# Write the new protein child class here
class protein(seq):
    def __init__(self, name, organism, sequence, type, size):
        super().__init__(name, organism, sequence, type)
        self.size = size

    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + str(self.size)
            + "\n"
            + self.sequence
            + "\n"
        )
        f.close()

    def mol_weight(self):
        total = 0
        for aa in self.sequence:
            total += aa_mol_weights.get(aa)
        return total


# # E. Using super(), make a child class of seq called nucleotide
# ## Make a new method called gc_content that calculates the percent of letters that are G or C and then prints the gc content percentage
#
# ## for the gc_content method, there are multiple ways to do it, but you may find this page helpful https://www.geeksforgeeks.org/python-count-occurrences-of-a-character-in-string/

# In[142]:


# Write the new nucleotide class here
class nucleotide(seq):

    def gc_content(self):
        total = len(self.sequence)
        g_count = self.sequence.count("G")
        c_count = self.sequence.count("C")
        percentage = ((g_count + c_count) / total) * 100

        print(percentage)


# # G. Using super(),make "Grandchild classes" DNA and RNA, which will be child classes of nucleotide. So in this analogy they would be like "grandchildren" of sequence class. Look at the homework instructions pdf for a schematic.

# ## G.1 For the DNA class, add a method called transcribe to transcribe the DNA to RNA and print the transcribed sequence (aka replace the Ts in the DNA sequence with Us)

# In[144]:


# Write the DNA class here
class DNA(nucleotide):

    # For the transcribe method
    # This site may be helpful https://www.geeksforgeeks.org/python-string-replace/

    def transcribe(self):
        transcribed = self.sequence.replace("T", "U")
        print(transcribed)

    def six_frames(self):
        frames = []

        frames.append(self.sequence)
        frames.append(self.sequence[1:])
        frames.append(self.sequence[2:])

        reverse_complement = self.reverse_complement()
        frames.append(reverse_complement)
        frames.append(reverse_complement[1:])
        frames.append(reverse_complement[2:])

        return frames

    def reverse_complement(self):
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        reverse = self.sequence[::-1]
        rc = ""
        for base in reverse:
            com = complement.get(base)
            rc += com
        return rc


# ## G.3 For the RNA class, add a method called start that will print the index of the start codon (AUG) in a RNA sequence

# In[146]:


# Write the RNA class here
class RNA(nucleotide):

    # For the start method
    # This site may be https://www.geeksforgeeks.org/python-string-find/
    def start(self):
        index = self.sequence.find("AUG")

    def translate(self):

        index = self.sequence.find("AUG")

        codons = []
        for i in range(index, len(self.sequence), 3):
            codon = self.sequence[i : i + 3]
            codons.append(codon)

        sequence = ""
        for codon in codons:
            aa = standard_code.get(codon)
            if aa == "*":
                break
            sequence += aa

        return sequence


# # Part B: test new methods

# In[148]:


uidA = DNA(
    name="uidA",
    organism="bacteria",
    sequence="CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA",
    type="DNA",
)

uidA.fasta_out()


# In[149]:

print(uidA.reverse_complement())
print(uidA.six_frames())


# In[150]:


print(uidA.transcribe())


uidA_RNA = RNA(
    name="uidA_RNA",
    organism="bacteria",
    sequence="CGCAUGUUACGUCCUGUAGAAACCCCAACCCGUGAAAUCAAAAAA",
    type="RNA",
)


# In[151]:


uidA_RNA.fasta_out()


# In[152]:


print(uidA_RNA.translate())


# In[153]:


uidA_protein = protein(
    name="uidA_protein",
    sequence="MLRPVETPTREIKK",
    organism="bacteria",
    type="protein",
    size="38",
)


# In[154]:


uidA_protein.fasta_out()


# In[155]:


print(uidA_protein.mol_weight())


# In[ ]:
