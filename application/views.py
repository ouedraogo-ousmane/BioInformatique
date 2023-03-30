from django.shortcuts import render
from django.http import HttpResponse
import numpy as np
from Bio.Seq import Seq
from Bio import pairwise2,SeqIO
from Bio.Align.Applications import ClustalwCommandline
import os
from django.conf import settings
# Create your views here.
def majuscule_view(request):
    if request.method == 'POST':
        chaine = request.POST.get('chaine', '')
        sequence = request.POST.get('sequence', '')
        seq1 = Seq(chaine)
        seq2 = Seq(sequence)
        alignement = pairwise2.align.globalxx(seq1,seq2)
        return render(request, 'majuscule.html', {'chaine_maj': alignement, 'chaine': seq2.complement()})
    else:
        return render(request, 'form.html')

    """Cette methode permet de determiner l'alignement globale ==
    de deux sequences seq1 et seq2
    """
def global_view(request):
    if request.method == 'POST':
        chaine = request.POST.get('chaine', '')
        sequence = request.POST.get('sequence', '')
        seq1 = Seq(chaine)
        seq2 = Seq(sequence)
        alignement = pairwise2.align.globalxx(seq1,seq2)
        return render(request, 'global.html', {'seq1': seq1, 'seq2': seq2,'resultat':alignement})
    else:
        return render(request, 'form.html')

def local_view(request):
    if request.method == 'POST':
        chaine = request.POST.get('chaine', '')
        sequence = request.POST.get('sequence', '')
        seq1 = Seq(chaine)
        seq2 = Seq(sequence)
        alignement = pairwise2.align.localxx(seq1,seq2)
        return render(request, 'local.html', {'seq1': seq1, 'seq2': seq2,'resultat':alignement})
    else:
        return render(request, 'form.html')


def global_alignGap_view(request):
    gap_open=-2
    gap_closed = -0.5
    match_score=1
    mismatch_penalty=-1
    if request.method == 'POST':
        chaine = request.POST.get('chaine', '')
        sequence = request.POST.get('sequence', '')
        seq1 = Seq(chaine)
        seq2 = Seq(sequence)
        alignement = pairwise2.align.globalms(seq1,seq2,match_score,mismatch_penalty,gap_open,gap_closed)
        return render(request, 'globalGap.html', {'seq1': seq1, 'seq2': seq2,'resultat':alignement})
    else:
        return render(request, 'form.html')


def clustal_Align_view(request):
    if request.method == 'POST':
        fichier = request.POST.get('fichier', '')
        print(fichier)
        chemin_rel = os.path.join(settings.MEDIA_ROOT, fichier)
        chemin_abs = os.path.abspath(chemin_rel)
        sequences = SeqIO.parse(chemin_abs, "fasta")
        clustalw_cline = ClustalwCommandline("clustalw", infile=chemin_abs)
        #stdout, stderr = clustalw_cline()
        alignement = AlignIO.read("fichier.aln", "clustal")
        return render(request, 'clustal.html', {'sequences':sequences,'resultat':alignement})
    else:
        return render(request, 'form1.html')

def arbre_WPGMA(request):
    pass

def needleman_wunsch_view(request):

    gap_penalty=-1
    match_score=1
    mismatch_penalty=-1
    if request.method == 'POST':
        chaine = request.POST.get('chaine', '')
        sequence = request.POST.get('sequence', '')
        s1 = Seq(chaine)
        s2 = Seq(sequence)
        # Création de la matrice de scores
        m, n = len(s1), len(s2)
        score_matrix = np.zeros((m+1, n+1))
        
        for i in range(1, m+1):
            score_matrix[i][0] = score_matrix[i-1][0] + gap_penalty
            
        for j in range(1, n+1):
            score_matrix[0][j] = score_matrix[0][j-1] + gap_penalty
            
        for i in range(1, m+1):
            for j in range(1, n+1):
                if s1[i-1] == s2[j-1]:
                    score_matrix[i][j] = score_matrix[i-1][j-1] + match_score
                else:
                    score_matrix[i][j] = max(score_matrix[i-1][j] + gap_penalty,
                                            score_matrix[i][j-1] + gap_penalty,
                                            score_matrix[i-1][j-1] + mismatch_penalty)
        
        # Récupération de l'alignement optimal
        align1, align2 = '', ''
        i, j = m, n
        
        while i > 0 and j > 0:
            if score_matrix[i][j] == score_matrix[i-1][j-1] + (match_score if s1[i-1] == s2[j-1] else mismatch_penalty):
                align1 = s1[i-1] + align1
                align2 = s2[j-1] + align2
                i -= 1
                j -= 1
            elif score_matrix[i][j] == score_matrix[i-1][j] + gap_penalty:
                align1 = s1[i-1] + align1
                align2 = '-' + align2
                i -= 1
            else:
                align1 = '-' + align1
                align2 = s2[j-1] + align2
                j -= 1
        
        while i > 0:
            align1 = s1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        
        while j > 0:
            align1 = '-' + align1
            align2 = s2[j-1] + align2
            j -= 1
        
        # return align1, align2
        return render(request, 'exact.html', {'seq1': s1, 'seq2': s1,'alignement1':align1,'alignement2':align2})
    else:
        return render(request, 'form.html')



