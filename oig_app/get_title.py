def make_title(mol_name, calculation_type, method_basis):

    calculation_type = calculation_type.replace('!', '')

    method_basis = method_basis.replace('!', '')
    title = str(mol_name) + str(calculation_type) + str(method_basis)
    title_string = f"#\n# {title}\n#"

    return title_string

    