from django.shortcuts import render
from django.http import HttpResponse
import numpy as np
from Bio.Seq import Seq
from Bio import pairwise2, SeqIO
from Bio.Align.Applications import ClustalwCommandline
import os
from django.conf import settings

from django.core.files.storage import FileSystemStorage
from Bio import Phylo

from django import forms
from io import StringIO
import re


class DndForm(forms.Form):
    file = forms.FileField()
# Create your views here.


def home_view(request):
    return render(request, 'index.html')


def dnd_view(request):
    if request.method == 'POST':
        form = DndForm(request.POST, request.FILES)
        if form.is_valid():
            # Récupérer le fichier DND soumis
            file = request.FILES['file']
            # Lire le contenu du fichier
            content = file.read().decode('utf-8')
            # Traiter le contenu pour générer l'arbre WPGMA
            tree = Phylo.read(StringIO(content), 'newick')
            # Rendre l'arbre WPGMA à l'aide d'un template
            return render(request, 'tree.html', {'tree': tree})
    else:
        form = DndForm()
    return render(request, 'upload1.html', {'form': form})


def clustal_Align_view(request):
    if request.method == 'POST' and request.FILES['fasta']:
        fasta_file = request.FILES['fasta']
        fs = FileSystemStorage()
        filename = fs.save(fasta_file.name, fasta_file)
        file_path = fs.path(filename)

        # Appel de Clustal pour effectuer l'alignement
        cline = ClustalwCommandline("ClustalW2", infile=file_path)
        stdout, stderr = cline()

        with open(file_path.replace('.fasta', '.aln'), 'r') as aln_file:
            aln_content = aln_file.read()

        return render(request, 'clustal.html', {'alignment': aln_content})

    return render(request, 'upload.html')

    """Cette methode permet de determiner l'alignement globale ==
    de deux sequences seq1 et seq2
    """


# Vérifie si la chaîne passée en paramètre ne contient que des nucléotides d'ADN
def is_valid_dna(seq):
    return bool(re.match('^[ATCG]*$', seq))

# Demande à l'utilisateur de saisir les deux séquences
# while True:
#     seq1 = input("Veuillez saisir la première séquence d'ADN : ")
#     if is_valid_dna(seq1):
#         break
#     print("La séquence saisie contient des caractères non-valides. Veuillez ne saisir que des nucléotides d'ADN (A, T, C ou G)")

# while True:
#     seq2 = input("Veuillez saisir la deuxième séquence d'ADN : ")
#     if is_valid_dna(seq2):
#         break
#     print("La séquence saisie contient des caractères non-valides. Veuillez ne saisir que des nucléotides d'ADN (A, T, C ou G)")

# # Effectue l'alignement local des deux séquences
# align1, align2 = smith_waterman(seq1, seq2)

# # Affiche les séquences alignées
# print(align1)
# print(align2)


def global_view(request):
    if request.method == 'POST':
        sequence_1 = request.POST.get('sequence_1', '')
        sequence_1 = sequence_1.strip().upper()

        sequence_2 = request.POST.get('sequence_2', '')
        sequence_2 = sequence_2.strip().upper()
        if(is_valid_dna(sequence_1) and is_valid_dna(sequence_2)):

            seq1 = Seq(sequence_1)
            seq2 = Seq(sequence_2)
            alignement = pairwise2.align.globalxx(seq1, seq2)
            return render(request, 'global-align.html', {'seq1': seq1, 'seq2': seq2, 'resultat': alignement})
        return render(request, 'global-align.html', {"error_check": True})
    else:
        return render(request, 'global-align.html')


def local_view(request):
    if request.method == 'POST':
        sequence_1 = request.POST.get('sequence_1', '')
        sequence_1 = sequence_1.strip().upper()

        sequence_2 = request.POST.get('sequence_2', '')
        sequence_2 = sequence_2.strip().upper()
        if(is_valid_dna(sequence_1) and is_valid_dna(sequence_2)):
            seq1 = Seq(sequence_1)
            seq2 = Seq(sequence_2)
            alignement = pairwise2.align.localxx(seq1, seq2)
            return render(request, 'global-align.html', {'seq1': seq1, 'seq2': seq2, 'resultat': alignement})
        return render(request, 'global-align.html', {"error_check": True})
    else:
        return render(request, 'local-align.html', {'post': False})


def global_alignGap_Affine_view(request):
    """Cette vue permet d'afficher l'alignement global avec des modeles de
       gap affine
    Args:
        request : la requete

    Returns:
        _type_: _description_
    """
    # gap_open = -2
    # gap_closed = -0.5
    # match_score = 1
    # mismatch_penalty = -1
    if request.method == 'POST':
        sequence_1 = request.POST.get('sequence_1', '')
        sequence_1 = sequence_1.strip().upper()

        sequence_2 = request.POST.get('sequence_2', '')
        sequence_2 = sequence_2.strip().upper()

        match_score = request.POST.get('match_score', '')
        match_score = int(match_score)

        mismatch_penalty = request.POST.get('mismatch_penalty', '')
        mismatch_penalty = int(mismatch_penalty)

        gap_open = request.POST.get('gap_open', '')
        gap_open = float(gap_open)

        gap_closed = request.POST.get('gap_closed', '')
        gap_closed = float(gap_closed)
        print(gap_closed, gap_open)

        seq1 = Seq(sequence_1)
        seq2 = Seq(sequence_2)

        if(is_valid_dna(sequence_1) and is_valid_dna(sequence_2)):
            seq1 = Seq(sequence_1)
            seq2 = Seq(sequence_2)
            alignement = pairwise2.align.globalms(
                seq1, seq2, match_score, mismatch_penalty, gap_open, gap_closed)
            return render(request, 'affine.html', {'post': False,
                                                   'match': match_score, 'missmatch': mismatch_penalty,
                                                   'gap_open': gap_open, 'gap_closed': gap_closed,
                                                   'seq1': seq1, 'seq2': seq2, 'resultat': alignement})
        return render(request, 'affine.html', {"error_check": True})
    else:
        return render(request, 'affine.html', {'post': False})


def global_alignGap_lineaire_view(request):
    """_summary_

    Args:
        request (_type_): _description_

    Returns:
        _type_: _description_
    """
    # match_score = 1
    # mismatch_penalty = -1
    if request.method == 'POST':
        sequence_1 = request.POST.get('sequence_1', '')
        sequence_1 = sequence_1.strip().upper()
        sequence_2 = request.POST.get('sequence_2', '')
        sequence_2 = sequence_2.strip().upper()
        match_score = request.POST.get('match_score', '')
        match_score = int(match_score)

        mismatch_penalty = request.POST.get('mismatch_penalty', '')
        mismatch_penalty = int(mismatch_penalty)
        if(is_valid_dna(sequence_1) and is_valid_dna(sequence_2)):
            seq1 = Seq(sequence_1)
            seq2 = Seq(sequence_2)
            alignement = pairwise2.align.globalmx(
                seq1, seq2, match_score, mismatch_penalty)
            return render(request, 'linear.html', {'post': False,
                                                   'match': match_score, 'missmatch': mismatch_penalty,
                                                   'seq1': seq1, 'seq2': seq2, 'resultat': alignement})
        return render(request, 'linear.html', {"error_check": True})
    else:
        return render(request, 'linear.html', {'post': False})


def multiple_align(request):
    seq = []
    if request.method == 'POST':
        if request.POST.get('subject') == 'Ajouter':
            sequence_1 = request.POST.get('sequence_1', '')
            seq.append(sequence_1)
            return render(request, 'test.html', {'resultat': seq})
        if request.POST.get('subject') == 'Convertir':
            return render(request, 'test.html', {'resultat': "Hello dear family"})
    return render(request, 'formulaire.html')


# def clustal_Align_view(request):
#     if request.method == 'POST':
#         fichier = request.POST.get('fichier', '')
#         print(fichier)
#         chemin_rel = os.path.join(settings.MEDIA_ROOT, fichier)
#         chemin_abs = os.path.abspath(chemin_rel)
#         sequences = SeqIO.parse(chemin_abs, "fasta")
#         clustalw_cline = ClustalwCommandline("clustalw", infile=chemin_abs)
#         #stdout, stderr = clustalw_cline()
#         alignement = AlignIO.read("fichier.aln", "clustal")
#         return render(request, 'clustal.html', {'sequences': sequences, 'resultat': alignement})
#     else:
#         return render(request, 'form1.html')


def arbre_WPGMA(request):

    if request.method == 'POST':
        fichier = request.POST.get('fichier', '')
        chemin_rel = os.path.join(settings.MEDIA_ROOT, fichier)
        chemin_abs = os.path.abspath(chemin_rel)
        tree = Phylo.read(chemin_abs, "newick")

        return render(request, 'clustal.html', {'resultat': Phylo.draw_ascii(tree)})
    else:
        return render(request, 'form1.html')


def needleman_wunsch_view(request):

    gap_penalty = -1
    match_score = 1
    mismatch_penalty = -1
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
                                             score_matrix[i][j-1] +
                                             gap_penalty,
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
        return render(request, 'exact.html', {'seq1': s1, 'seq2': s1, 'alignement1': align1, 'alignement2': align2})
    else:
        return render(request, 'form.html')
