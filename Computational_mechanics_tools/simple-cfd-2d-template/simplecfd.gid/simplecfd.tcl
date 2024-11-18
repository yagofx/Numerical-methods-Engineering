
#define the procedures in a separated namespace named SimpleCFD
namespace eval SimpleCFD {
    variable conditions_list
}

proc SimpleCFD::Init { } {
    variable element_id
    variable conditions_list
    variable boundary_conditions_list

    # TODO: Add here the id of your fluid from the spd.
    # if your spd has: <condition n="fluid" ... > then you should set element_id "fluid"
    set element_id "YOUR FLUID CONDITION ID"
    # TODO: Add here the list of ids of tyour boundary conditions from the spd
    # if your spd has: <condition n="inlet" ... > then you should add "inlet" to the list
    # if your spd has: <condition n="outlet" ... > then you should add "outlet" to the list
    # ...
    set boundary_conditions_list [list "YOUR_BOUNDARY_CONDITION_ID_1" "YOUR_BOUNDARY_CONDITION_ID_2"]

    set conditions_list [concat $element_id $boundary_conditions_list]
}

# DO NOT TOUCH ANYTHING BELOW THIS LINE
proc AfterWriteCalcFileGIDProject { filename errorflag } {
    if { ![info exists gid_groups_conds::doc] } {
        WarnWin [= "Error: data not OK"]
        return
    }
    set err [SimpleCFD::ValidateData]
    set ret ""
    if { $err > 0} {
        GiD_WriteCalculationFile end
        WarnWin [= "Error when preparing data for analysis (%s)" $::errorInfo]
        set ret -cancel-
    }
    if { $ret == "-cancel-" } {
        return $ret
    }
    set err [catch { SimpleCFD::WriteModelFile $filename } ret]
    if { $err } {
        GiD_WriteCalculationFile end
        WarnWin [= "Error when preparing data for analysis (%s)" $::errorInfo]
        set ret -cancel-
    }
    if { $ret == "-cancel-" } {
        return $ret
    }
    set fname [file join [file dirname $filename] "ProjectParameters.json"]
    set err [catch { SimpleCFD::WriteParametersFile $fname} ret]
    if { $err } {
        GiD_WriteCalculationFile end
        WarnWin [= "Error when preparing data for analysis (%s)" $::errorInfo]
        set ret -cancel-
    }
    if { $ret == "-cancel-" } {
        return $ret
    }
    set fname [file join [file dirname $filename] "FluidMaterials.json"]
    set err [catch { SimpleCFD::WriteMaterialsFile $fname} ret]
    if { $err } {
        GiD_WriteCalculationFile end
        WarnWin [= "Error when preparing data for analysis (%s)" $::errorInfo]
        set ret -cancel-
    }

    return $ret
}

proc SimpleCFD::ValidateData {} {
    set err 0
    # Check that the model is meshed
    if { [GiD_Info Mesh] eq 0 } {
        WarnWin [= "Error: model not meshed"]
        incr err
    }
    # Check that only 1 fluid is defined
    set xpath "//condition\[@n='$SimpleCFD::element_id'\]/group"
    set xml_nodes [[$::gid_groups_conds::doc documentElement] selectNodes $xpath]
    if {[llength $xml_nodes] ne 1} {
        WarnWin "Error: only one fluid can be defined"
        incr err
    }
    # Check that the fluid is meshed
    set fluid_id [[lindex $xml_nodes 0] @n]
    if { [GiD_EntitiesGroups get $fluid_id elements -count] eq 0 } {
        WarnWin [= "Error: fluid not meshed. remember to mesh the fluid"]
        incr err
    }

    # Check that the lines are meshed
    foreach condition $SimpleCFD::boundary_conditions_list {
        set xpath "//condition\[@n='$condition'\]/group"
        set xml_nodes [[$::gid_groups_conds::doc documentElement] selectNodes $xpath]
        
        foreach group $xml_nodes {
            set group_id [$group @n]
            if { [GiD_EntitiesGroups get $group_id elements -count] eq 0 } {
                WarnWin [= "Error: group $group_id not meshed. Remember to mesh the lines."]
                incr err
            }
        }
    }

    return $err
}

proc InitGIDProject { dir } {
    
    # Load my scripts
    source [file join $dir "write.tcl"]
    source [file join $dir "write_parameters.tcl"]
    source [file join $dir "write_materials.tcl"]
    
    # set environment variables for python
    set python_path [GiD_Python_GetPythonExe]
    set ::env(python_path) $python_path
    set ::env(python_home) [file dirname $python_path]

    gid_groups_conds::open_conditions menu
    SimpleCFD::Init
}
