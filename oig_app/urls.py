from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('inputs/', views.inputs, name='inputs'),
    path('jmol_view', views.jmol_view, name='jmol_view'),
    path('submitted/', views.submitted, name='submitted'),
]   