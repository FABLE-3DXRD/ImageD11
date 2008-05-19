if {![package vsatisfies [package provide Tcl] 8]} {return}
package ifneeded Togl 1.6.0 \
    [list load [file join $dir Togl.dll] Togl]