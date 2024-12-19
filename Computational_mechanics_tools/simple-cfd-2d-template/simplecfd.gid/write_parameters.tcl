proc SimpleCFD::WriteParametersFile {fname} {

    # Load templates
    set modeler_condition_template [SimpleCFD::ReadTemplateFile modeler_condition]
    set modeler_element_template [SimpleCFD::ReadTemplateFile modeler_element]
    set inlet_template [SimpleCFD::ReadTemplateFile inlet]
    set outlet_template [SimpleCFD::ReadTemplateFile outlet]
    set noslip_template [SimpleCFD::ReadTemplateFile noslip]
    set slip_template [SimpleCFD::ReadTemplateFile slip]
    set main_template [SimpleCFD::ReadTemplateFile main]

    GiD_WriteCalculationFile init $fname

    # remove the folders in path
    set case_name [file tail [GiD_Info Project ModelName]]

    # Time parameters
    set start_time 0.0
    # Get the end time from the spd file
    set end_time_xpath "//value\[@n='end_time'\]"
    set end_time [get_domnode_attribute [[$::gid_groups_conds::doc documentElement] selectNodes $end_time_xpath] v]
    
    # TODO: Do the same for the delta time
    set time_step_xpath "//value\[@n='delta_time'\]"
    set time_step [get_domnode_attribute [[$::gid_groups_conds::doc documentElement] selectNodes $time_step_xpath] v]
    # Check the variable name in the main template

    #Strategy parameters
    set max_iterations_xpath "//value\[@n='max_iterations'\]"
    set max_iterations [get_domnode_attribute [[$::gid_groups_conds::doc documentElement] selectNodes $max_iterations_xpath] v]

    set relative_velocity_tolerance_xpath "//value\[@n='rel_vel_tol'\]"
    set relative_velocity_tolerance [get_domnode_attribute [[$::gid_groups_conds::doc documentElement] selectNodes $relative_velocity_tolerance_xpath] v]

    set absolute_velocity_tolerance_xpath "//value\[@n='abs_vel_tol'\]"
    set absolute_velocity_tolerance [get_domnode_attribute [[$::gid_groups_conds::doc documentElement] selectNodes $absolute_velocity_tolerance_xpath] v]

    set relative_pressure_tolerance_xpath "//value\[@n='rel_press_tol'\]"
    set relative_pressure_tolerance [get_domnode_attribute [[$::gid_groups_conds::doc documentElement] selectNodes $relative_pressure_tolerance_xpath] v]

    set absolute_pressure_tolerance_xpath "//value\[@n='abs_press_tol'\]"
    set absolute_pressure_tolerance [get_domnode_attribute [[$::gid_groups_conds::doc documentElement] selectNodes $absolute_pressure_tolerance_xpath] v]



    


    # Conditions
    set list_skin_parts ""
    foreach group [SimpleCFD::GetSkinModelpartsGroups] {
        set group_name [SimpleCFD::TransformGroupName $group]
        append list_skin_parts "\"$group_name\","
    }
    if { $list_skin_parts ne "" } {
        set list_skin_parts [string trim $list_skin_parts ","]
    }

    set group [SimpleCFD::GetBodyModelpartGroup] 
    set group_name [SimpleCFD::TransformGroupName $group]
    set volume_model_part_name $group_name
    
    # Fill this list with the processes that are going to be used in the boundary conditions
    # TODO: Add the other boundary conditions
    # I provide the outlet as example, add the other ones checking in the templates if there are any parameters to be filled
    # You can find examples on how to get and fill the parameters in the ent_time above, or in the write_materials.tcl script
    set boundary_conditions_process_list ""
    
    # TODO: Add the inlet
    foreach group [SimpleCFD::GetUsedGroups "automatic_inlet_velocity"] {
        set inlet_value_xpath "//condition\[@n='automatic_inlet_velocity'\]/group/value\[@n='velocity'\]"
        set inlet_value [get_domnode_attribute [[$::gid_groups_conds::doc documentElement] selectNodes $inlet_value_xpath] v]
        set modelpart_name [SimpleCFD::TransformGroupName $group]
        set result [subst -nobackslashes -nocommands $inlet_template]
        append boundary_conditions_process_list "$result,"
    }

    # Outlets
    foreach group [SimpleCFD::GetUsedGroups "outlet_pressure"] {
        set outlet_value_xpath "//condition\[@n='outlet_pressure'\]/group/value\[@n='pressure'\]"
        set outlet_value [get_domnode_attribute [[$::gid_groups_conds::doc documentElement] selectNodes $outlet_value_xpath] v]
        set modelpart_name [SimpleCFD::TransformGroupName $group]
        set result [subst -nobackslashes -nocommands $outlet_template]
        append boundary_conditions_process_list "$result,"
    }
    # TODO: Add the slip condition
    foreach group [SimpleCFD::GetUsedGroups "slip"] {
        set modelpart_name [SimpleCFD::TransformGroupName $group]
        set result [subst -nobackslashes -nocommands $slip_template]
        append boundary_conditions_process_list "$result,"
    }
    # TODO: Add the noslip condition
    foreach group [SimpleCFD::GetUsedGroups "no_slip"] {
        set modelpart_name [SimpleCFD::TransformGroupName $group]
        set result [subst -nobackslashes -nocommands $noslip_template]
        append boundary_conditions_process_list "$result,"
    }
    
    # remove the last comma
    if { $boundary_conditions_process_list ne "" } {
        set boundary_conditions_process_list [string trim $boundary_conditions_process_list ","]
    }

    # Modelers section
    set modeler_elements_list ""
    set group [SimpleCFD::GetBodyModelpartGroup] 
    set modelpart_name [SimpleCFD::TransformGroupName $group]
    set element_name "Element2D3N"
    set result [subst -nobackslashes -nocommands $modeler_element_template]
    append modeler_elements_list "$result,"
    if { $modeler_elements_list ne "" } {
        set modeler_elements_list [string trim $modeler_elements_list ","]
    }

    set modeler_conditions_list ""
    foreach group [SimpleCFD::GetSkinModelpartsGroups] {
        set modelpart_name [SimpleCFD::TransformGroupName $group]
        set condition_name "WallCondition2D2N"
        set result [subst -nobackslashes -nocommands $modeler_condition_template]
        append modeler_conditions_list "$result,"
    }
    if { $modeler_conditions_list ne "" } {
        set modeler_conditions_list [string trim $modeler_conditions_list ","]
    }

    set result [subst -nobackslashes -nocommands $main_template]

    GiD_WriteCalculationFile puts $result
    GiD_WriteCalculationFile end
}


# Function to get the list of names of the groups used in the boundary conditions
proc SimpleCFD::GetSkinModelpartsGroups { } {
    set groups [list ]

    foreach condition $SimpleCFD::boundary_conditions_list {
        set groups [concat $groups [SimpleCFD::GetUsedGroups $condition]]
    }

    set all [lsort -unique $groups]
    return $all
}

# Function to get the name of used in the fluid condition
proc SimpleCFD::GetBodyModelpartGroup { } {
    return [lindex [SimpleCFD::GetUsedGroups $SimpleCFD::element_id] 0]
}

proc SimpleCFD::ReadTemplateFile { name } {
    set template_file [file join [GiD_Info problemtypepath] templates "$name.json"]
    set fp [open $template_file r]
    set template [read $fp]
    close $fp
    return $template
}

