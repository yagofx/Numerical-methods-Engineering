proc SimpleCFD::WriteMaterialsFile {fname} {
    GiD_WriteCalculationFile init $fname
    
    # Read the material template
    set materials_template [SimpleCFD::ReadTemplateFile material]

    # Get the density from the spd file
    set density_xpath "//condition\[@n='fluid'\]/group/value\[@n='density'\]"
    set density [get_domnode_attribute [[$::gid_groups_conds::doc documentElement] selectNodes $density_xpath] v]
    
    # TODO: Do the same for the dynamic viscosity
    # Check the variable name in the materials template
    set dynamic_viscosity_xpath "//condition\[@n='fluid'\]/group/value\[@n='viscosity'\]"
    set dynamic_viscosity [get_domnode_attribute [[$::gid_groups_conds::doc documentElement] selectNodes $dynamic_viscosity_xpath] v]

    set group [SimpleCFD::GetBodyModelpartGroup]

    set element_type [SimpleCFD::GetElementTypeUsedByGroup $group]
    if { $element_type == "Triangle"} {
        set NewtonianLaw "Newtonian2DLaw"
    }
    if { $element_type == "Tetrahedra"} {
        set NewtonianLaw "Newtonian3DLaw"
    }

    set group_name [SimpleCFD::TransformGroupName $group]
    set model_part_name $group_name
    GiD_WriteCalculationFile puts [subst -nobackslashes -nocommands $materials_template]
    GiD_WriteCalculationFile end
}

proc SimpleCFD::GetElementTypeUsedByGroup {group} {
    return [lindex [GiD_Mesh get element [objarray get [GiD_EntitiesGroups get $group elements] 0] ] 1]
}   
