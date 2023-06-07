from django.urls import path
from . import views

print('Inside urls.py')

urlpatterns = [
    path('', views.index, name='index'),
    path('inputs/', views.inputs, name='inputs'),
    path('mol_view/', views.mol_view, name='mol_view'),
    path('submitted/', views.submitted, name='submitted'),
]   