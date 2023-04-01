from django.urls import path
from django.conf import settings
from django.conf.urls.static import static
from .views import (global_view,local_view,
needleman_wunsch_view,multiple_align,global_alignGap_Affine_view,clustal_Align_view,global_alignGap_lineaire_view,dnd_view,upload_file,arbre_WPGMA)

urlpatterns = [
    path('global/', global_view, name='global'),
    path('local/', local_view, name='local'),
    path('multiple/', multiple_align, name='multiple'),
    path('exacte/', needleman_wunsch_view, name='exacte'),
    path('globalGapAffine/', global_alignGap_Affine_view, name='globalGapAffine'),
    path('globalGapLineaire/', global_alignGap_lineaire_view, name='globalGapLineaire'),
    path('clustal/', upload_file, name='clustal'),
    path('arbre/', dnd_view, name='arbre'),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
