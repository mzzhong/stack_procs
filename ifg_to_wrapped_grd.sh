#!/bin/bash
#filename=diff_20131222_20131223_20131019_20131020
#filename=diff_20131222_20131223_20130901_20130902
filename=diff_20140814_20140815_20140729_20140730

gdal_translate -of GMT -b 1 ${filename}_m.off ${filename}_m.grd
grdmath ${filename}_m.grd 0 NAN = ${filename}_m_masked.grd 
grdmath ${filename}_m_masked.grd 0.02 MOD = ${filename}_m_masked_wrapped.grd
grdmath ${filename}_m_masked_wrapped.grd 0.01 SUB = ${filename}_m_masked_wrapped.grd

#grdmath ${filename}_m_wrapped.grd 0 NAN = ${filename}_m_wrapped_masked.grd 
