DESIGN=/proj/ipeng2-nobkup/sprusty/proj/vcu/isp_integration/vivado/xilinx_isp_explore_AltSpreadLogic_medium_2019.1_withrevsionlikeinterface/zcu104_vcu_broadcast
PROJ=`find $DESIGN -name *.xpr -printf "%f\n"`
echo ${PROJ}

USERNAME=${PROJ%.xpr*}
echo $USERNAME


