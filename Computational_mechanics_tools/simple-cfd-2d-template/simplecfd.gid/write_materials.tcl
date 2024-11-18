proc SimpleCFD::WriteMaterialsFile {fname} {
    GiD_WriteCalculationFile init $fname
    
    # Read the material template
    set materials_template [SimpleCFD::ReadTemplateFile material]

    # Get the density from the spd file
    set density_xpath "//condition\[@n='fluid'\]/group/value\[@n='density'\]"
    set density [get_domnode_attribute [[$::gid_groups_conds::doc documentElement] selectNodes $density_xpath] v]
    
    # TODO: Do the same for the dynamic viscosity
    # Check the variable name in the materials template


    

    set group [SimpleCFD::GetBodyModelpartGroup] 
    set group_name [SimpleCFD::TransformGroupName $group]
    set model_part_name $group_name
    GiD_WriteCalculationFile puts [subst -nobackslashes -nocommands $materials_template]
    GiD_WriteCalculationFile end
}
