document {
	Key => {(action, RingOfInvariants)},
	
	Headline => "the group action that produced a ring of invariants",
	
	Usage => "action S",
	
	Inputs => {
    "S" => RingOfInvariants => {"of the group action being returned"},
  },
	
	Outputs => {
		GroupAction => {"the action that produced the ring of invariants in the input"}
	},

	PARA {
    "This example shows how to recover the action of a
    torus that produced a certain ring of invariants.
    Note other action types are possible and may produce
    a different looking output."
  },
    	
	EXAMPLE {
		"R = QQ[x_1..x_4]",
		"T = diagonalAction(matrix {{0,1,-1,1},{1,0,-1,-1}}, R)",
		"S = R^T",
		"action S"
  }
}