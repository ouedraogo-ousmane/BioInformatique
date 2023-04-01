from django.urls import path
from django.conf import settings
from django.conf.urls.static import static
from . import views

urlpatterns = [
    path('', views.home_view, name='home'),
    path('linear-gap/', views.global_view, name='linear'),
    path('affine-gap/', views.global_view, name='affine'),
    path('global/', views.global_view, name='global'),
    path('local/', views.local_view, name='local'),
    path('exacte/', views.needleman_wunsch_view, name='exacte'),
    path('globalGap/', views.global_alignGap_view, name='globalGap'),
    path('clustal/', views.upload_file, name='clustal'),
    path('arbre/', views.dnd_view, name='arbre'),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
