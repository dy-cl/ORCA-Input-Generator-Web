# views.py
from django.shortcuts import render, redirect
from django.utils.safestring import mark_safe

from .get_xyz import to_xyz
from .get_title import make_title
from .draw_mol import draw_mol

import base64
from io import BytesIO
import re

#HOME PAGE
##########################################
def index(request):

    return render(request, 'index.html')
##########################################

#INPUTS PAGE
##########################################
def inputs(request):
    if request.method == 'POST':
        calculation_type = request.POST.get('calculation_type')
        mol_input_format = request.POST.get('mol_input_format')
        method = request.POST.get('method')
        basis_set = request.POST.get('basis_set')
        mol_string = request.POST.get('mol_string')
        mol_name = request.POST.get('mol_name')
        coordinate_type = request.POST.get('coordinate_type')
        charge = request.POST.get('charge')
        spin = request.POST.get('spin')
        aux_basis_set = request.POST.get('aux_basis_set')
        excited_state_method = request.POST.get("excited_state_method")

        #Input validation and error checking
        ############################################# 
        empty_fields = [
            'mol_string',
            'charge',
            'spin',
            'mol_name',
        ]

        if any(not request.POST.get(field) for field in empty_fields):
            return render(request, 'inputs.html', {'empty_fields_error': 'Please fill in all sections of the form.'})

        if len(mol_string) > 1000:

            return render(request, 'inputs.html', {'mol_string_error': 'Molecular string cannot be more than 1000 characters.'})
        
        if len(mol_name) > 200:

            return render(request, 'inputs.html', {'mol_name_error': 'Molecule name cannot be more than 200 characters.'})
        
        if mol_input_format == 'smiles' and not re.match(r'^[^J][0-9BCOHNSOPrIFla@+\-\[\]\(\)\\\/%=#$,.~&!]+$', mol_string): 
            return render(request, 'inputs.html', {'invalid_smiles_error': 'Invalid SMILES string.'})
        
        if mol_input_format == 'inchi' and not re.match(r'^InChI\=1S?\/[A-Za-z0-9\.]+(\+[0-9]+)?(\/[cnpqbtmsih][A-Za-z0-9\-\+\(\)\,\/\?\;\.]+)*$', mol_string):
            return render(request, 'inputs.html', {'invalid_inchi_error': 'Invalid InChi string.'})
        
        if not charge.isdigit():
            return render(request, 'inputs.html', {'charge_error': 'Charge must be a numeric value.'})

        if not spin.isdigit():
            return render(request, 'inputs.html', {'spin_error': 'Spin must be a numeric value.'})
        
        ###############################################

        method_basis = '!' + ' ' + method + ' ' + basis_set

        xyz_content = to_xyz(mol_string, mol_input_format)
        title = make_title(mol_name, calculation_type, method_basis)
        csc = '*' + ' ' + str(coordinate_type) + ' ' + str(charge) + ' ' + str((int(spin)*2) + 1)

        request.session['xyz_content'] = xyz_content
        request.session['title'] = title
        request.session['method_basis'] = method_basis
        request.session['calculation_type'] = calculation_type
        request.session['csc'] = csc
        request.session['aux_basis_set'] = aux_basis_set
        request.session['excited_state_method'] = excited_state_method
        request.session['mol_string'] = mol_string
        request.session['mol_input_format'] = mol_input_format

        return redirect('mol_view')

    return render(request, 'inputs.html')
##########################################

#REVIEW 3D STRUCTURE REQUEST PAGE
##########################################
def mol_view(request):

    xyz_content = request.session.get('xyz_content')

    mol_image = draw_mol(xyz_content)

    mol_image = base64.b64encode(mol_image.getvalue()).decode('utf-8')

    context = {
        'mol_image': mol_image,
    }

    return render(request, 'mol_view.html', context)
##########################################

#GENERATED INPUT FILE PAGE
##########################################
def submitted(request):

    xyz_content = request.session.get('xyz_content')
    title = request.session.get('title')
    method_basis = request.session.get('method_basis')
    calculation_type = request.session.get('calculation_type')
    csc = request.session.get('csc')
    aux_basis_set = request.session.get('aux_basis_set')
    excited_state_method = "%" + request.session.get('excited_state_method').lower()

    context = {
        'xyz_content': xyz_content,
        'title': title,
        'method_basis': method_basis,
        'calculation_type': calculation_type,
        'csc': csc,
        'aux_basis_set': aux_basis_set,
        'excited_state_method': excited_state_method,
    }

    print(calculation_type)

    return render(request, 'submitted.html', context)
##########################################