from django.urls import path
from . import views
# (majuscule_view,global_view,local_view,
# needleman_wunsch_view,global_alignGap_view,clustal_Align_view)

urlpatterns = [
    path('', views.home_view, name='home'),
    path('majuscule/', views.majuscule_view, name='majuscule'),
    path('global/', views.global_view, name='global'),
    path('local/', views.local_view, name='local'),
    path('exacte/', views.needleman_wunsch_view, name='exacte'),
    path('globalGap/', views.global_alignGap_view, name='globalGap'),
    path('clustal/', views.clustal_Align_view, name='clustal'),
]
