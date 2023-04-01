from django.urls import path
from django.conf import settings
from django.conf.urls.static import static
from .views import (global_view,local_view,
needleman_wunsch_view,global_alignGap_view,clustal_Align_view,dnd_view,upload_file,arbre_WPGMA)

urlpatterns = [
    path('global/', global_view, name='global'),
    path('local/', local_view, name='local'),
    path('exacte/', needleman_wunsch_view, name='exacte'),
    path('globalGap/', global_alignGap_view, name='globalGap'),
    path('clustal/', upload_file, name='clustal'),
    path('arbre/', dnd_view, name='arbre'),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
