function trap_ctrl_c()
{
	echo "Enter"
	exit 2
}

trap trap_ctrl_c 2

sleep 10000

