from django.urls import path
from .views import (majuscule_view,global_view,local_view,
needleman_wunsch_view,global_alignGap_view,clustal_Align_view)

urlpatterns = [
    path('majuscule/', majuscule_view, name='majuscule'),
    path('global/', global_view, name='global'),
    path('local/', local_view, name='local'),
    path('exacte/', needleman_wunsch_view, name='exacte'),
    path('globalGap/', global_alignGap_view, name='globalGap'),
    path('clustal/', clustal_Align_view, name='clustal'),
]
