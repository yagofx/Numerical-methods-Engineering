
###################################################################################
#      print data in the .dat calculation file
################################################################################### 

proc SimpleCFD::WriteModelFile { filename } {
    GiD_WriteCalculationFile init $filename

    # write the header
    SimpleCFD::WriteCalculationFileHeader

    # write the nodes
    SimpleCFD::WriteNodalCoordinates

    # write the geometries
    SimpleCFD::WriteGeometries

    # write the modelparts
    SimpleCFD::WriteSubmodelparts

    GiD_WriteCalculationFile end

}

proc SimpleCFD::WriteCalculationFileHeader { } {
    # TODO: Write header




}

proc SimpleCFD::WriteNodalCoordinates { } {
    # TODO: Write nodes and coordinates




}

proc SimpleCFD::WriteGeometries { } {
    # foreach condition in the list of conditions defined in the simplecfd.tcl file
    foreach condition $SimpleCFD::conditions_list {
        # Write geometries for the condition
        SimpleCFD::WriteGeometriesOnCondition $condition
    }
}

# Write geometries for a given condition
proc SimpleCFD::WriteGeometriesOnCondition { condition_name } {
    # Get the list of groups used by the condition
    set used_groups [SimpleCFD::GetUsedGroups $condition_name]
    # For each group used by the condition
    foreach group $used_groups {
        # Get the element type used by the group (linear, triangle, etc.)
        set element_type [SimpleCFD::GetElementTypeUsedByGroup $group]
        # Get the element type that kratos expects (Line2D2, Triangle2D3, etc.)
        set kratos_element_type [SimpleCFD::TransformElementTypeName $element_type]
        # TODO: Write the geometry information (Check mdpa example -> Begin Geometries)




    }
}

proc SimpleCFD::WriteSubmodelparts { } {
    # foreach condition in the list of conditions defined in the simplecfd.tcl file
    foreach condition $SimpleCFD::conditions_list {
        # Write submodelparts for the condition
        SimpleCFD::WriteSubmodelpartsOnCondition $condition
    }
}


# Write submodelparts for a given condition
proc SimpleCFD::WriteSubmodelpartsOnCondition { condition_name } {
    # Get the list of groups used by the condition
    set used_groups [SimpleCFD::GetUsedGroups $condition_name]
    
    # For each group used by the condition
    foreach group $used_groups {
        # Get the group name transformed (spaces replaced by underscores)
        set group_name [SimpleCFD::TransformGroupName $group]
        
        # TODO: Write the submodelpart information (Check mdpa example -> Begin Submodelparts)





    }
}

# Do not modify the following functions

proc SimpleCFD::GetElementTypeUsedByGroup {group} {
    return [lindex [GiD_Mesh get element [objarray get [GiD_EntitiesGroups get $group elements] 0] ] 1]
}

# Get the list of groups used by a condition
proc SimpleCFD::GetUsedGroups { {condition_name ""} } {
    set used_groups [list ]
    set xpath "//group"
    if {$condition_name ne ""} {
        set xpath "//condition\[@n='$condition_name'\]/group"
    }
    set xml_nodes [[$::gid_groups_conds::doc documentElement] selectNodes $xpath]
    foreach group $xml_nodes {
        set group_name [get_domnode_attribute $group n]
        # Only write the groups that are used, and only once
        if { [lsearch $used_groups $group_name] != -1 } {
            continue
        }
        lappend used_groups $group_name
    }
    return $used_groups
}

proc SimpleCFD::TransformGroupName { group } {
    # Replace spaces with underscores
    set group [string map {" " "_"} $group]
    return $group
}

proc SimpleCFD::TransformElementTypeName { element_type } {
    set result ""
    if {$element_type eq "Linear"} {set result "Line2D2"}
    if {$element_type eq "Triangle"} {set result "Triangle2D3"}
    return $result
}
