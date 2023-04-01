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
from django.shortcuts import render
from io import StringIO


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


def upload_file(request):
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


def global_view(request):
    if request.method == 'POST':
        sequence_1 = request.POST.get('sequence_1', '')
        sequence_2 = request.POST.get('sequence_2', '')
        seq1 = Seq(sequence_1)
        seq2 = Seq(sequence_2)
        alignement = pairwise2.align.globalxx(seq1, seq2)
        return render(request, 'global-align.html', {'post': True, 'seq1': seq1, 'seq2': seq2, 'resultat': alignement})
    else:
        return render(request, 'global-align.html', {'post': False})


def local_view(request):
    if request.method == 'POST':
        sequence_1 = request.POST.get('sequence_1', '')
        sequence_2 = request.POST.get('sequence_2', '')
        seq1 = Seq(sequence_1)
        seq2 = Seq(sequence_2)
        alignement = pairwise2.align.localxx(seq1, seq2)
        print(alignement)
        return render(request, 'local-align.html', {'post': True, 'seq1': seq1, 'seq2': seq2, 'resultat': alignement})
    else:
        return render(request, 'local-align.html', {'post': False})


def global_alignGap_Affine_view(request):
    gap_open = -2
    gap_closed = -0.5
    match_score = 1
    mismatch_penalty = -1
    if request.method == 'POST':
        chaine = request.POST.get('chaine', '')
        sequence = request.POST.get('sequence', '')
        seq1 = Seq(chaine)
        seq2 = Seq(sequence)
        alignement = pairwise2.align.globalms(
            seq1, seq2, match_score, mismatch_penalty, gap_open, gap_closed)
        return render(request, 'globalGap.html', {'seq1': seq1, 'seq2': seq2, 'resultat': alignement})
    else:
        return render(request, 'form.html')


def global_alignGap_lineaire_view(request):
    gap_open = -2
    gap_closed = -0.5
    match_score = 1
    mismatch_penalty = -1
    if request.method == 'POST':
        chaine = request.POST.get('chaine', '')
        sequence = request.POST.get('sequence', '')
        seq1 = Seq(chaine)
        seq2 = Seq(sequence)
        alignement = pairwise2.align.globalmx(
            seq1, seq2, match_score, mismatch_penalty)
        return render(request, 'globalGap.html', {'seq1': seq1, 'seq2': seq2, 'resultat': alignement})
    else:
        return render(request, 'form.html')


def multiple_align(request):
    if request.method == 'POST':
        return render(request, 'test.html', {'resultat': request})
    return render(request, 'formulaire.html')


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
        return render(request, 'clustal.html', {'sequences': sequences, 'resultat': alignement})
    else:
        return render(request, 'form1.html')


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
