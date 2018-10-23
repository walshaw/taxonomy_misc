#!/usr/bin/perl

use strict;
use warnings;

use Carp;

use List::Util qw( first );

my $default_max_iterations = 100;

my $delimiter        = q{|};
#my $delimiter_qr     = qr{$delimiter}xms;
my $delimiter_quoted = quotemeta $delimiter;

#print qq{$delimiter_qr\n}; exit;

my $usage = qq{Usage:\n\n$0 TAXONOMY_NODES_FILE ANCESTRAL_NODE_ID [ MAX_ITERATIONS  [ REPORT_FREQUENCY ] ]

Probably want to use nodes.dmp from NCBI Taxonomy, as the first arg.
ANCESTRAL_NODE_ID is the node at the very top of the tree of interest.
The file is repeatedly read* until no further nodes are added to the clade,
or until the file has been read MAX_ITERATIONS times.

The point is that it cannot be assumed that the nodes are listed in a
'hierarchical order', i.e. it cannot be assumed that every node is listed
after its parent, grandparent and all ancestral nodes.

*Actually, although the file may be large - currently nodes.dmp is 1.3 million
lines long - it is still read only once (only the first two fields are read,
which both contain integer values), and then the data iterated over.
};

my $taxonomy_nodes_file = shift or croak $usage;

my $ancestral_node      = shift or croak $usage;

my $max_iterations      = shift || $default_max_iterations; # so this can't be 0

my $report_frequency    = shift || 20,000;

my $universal_clade_node = 1; # the ancestor of everything else;

my @parent_of; # indices correspond to node IDs; values are node IDs of parent

my @node_in_clade; # indices correspond to node IDs; values are boolean; if true then
                   # the node is part of the required clade; if defined but false
                   # the node is known to *not* be part of the clade; if undef then
                   # it's yet to be tested

open my $fh, '<', $taxonomy_nodes_file or croak qq{couldn't read nodes file '$taxonomy_nodes_file'};

my @clade_nodes = ($ancestral_node);

my $n_misordered_lines = 0;

LINE:
while (defined (my $line = <$fh>)) {

    my @fields = split m{ \s* $delimiter_quoted \s* }xms, $line;

    my ($node_id, $parent_id) = @fields[0,1];

    $parent_of[$node_id] = $parent_id;

=pod

    In most cases, a node's parent's entry will have been read *before* the
    (child) node itself's entry. But that is not guaranteed. Even if the
    nodes.dmp file has been sorted in ascending order of ordinal rank, some
    mis-ordering is possible because of all the 'no rank' entries, whose
    ordinal rank will have been merely inferred.

    The following tests for whether the current node's parent has been read
    yet; the way to do that is to check whether the parent node does itself
    have its parent node defined yet (reading a node's entry is one and the
    same as assigning its parent).

=cut

    if (!defined $parent_of[$parent_id]) {

        print STDERR qq{node $node_id has parent $parent_id; but entry for $parent_id has not been read yet\n};
        $n_misordered_lines++;
    }

#    if (first { $_ == $parent_id } @clade_nodes) {
#        push @clade_nodes, $node_id;
#        print STDERR qq{iteration 1: adding $node_id (parent is $parent_id)\n};
#    }
}

close $fh or croak qq{couldn't close nodes file '$taxonomy_nodes_file'};

print STDERR qq{$n_misordered_lines instances of nodes listed before their parent node\n};

NODE:
for my $node_id (0 .. $#parent_of) {

    next NODE if (!defined $parent_of[$node_id]); # some node_ids are unused in the nodes.dmp file

    my $parent_id = $parent_of[$node_id];

    my ($terminal_node, @lineage) = get_lineage($node_id);
    printf qq{node %8d\tIN CLADE:%s\tterminal node %8d; lineage: %s\n},
              $node_id, ($terminal_node == $ancestral_node) ? q{YES} : q{NO}, $terminal_node,
              join(q{, }, @lineage);
}




exit 0;

sub get_lineage {

=pod
    Starting with the specified node, an upwards lineage is assembled (parent; grandparent etc)
    until one of these four terminal conditions is met:

        the next node reached is already known to be a member of the clade of interest; or

        the externally-specified node of interest is reached (the ultimate ancestral
        node of the clade of interest); or

        the next node reached is already known to NOT be a member of the clade of interest; or

        the universal root node is reached.

    Either of the two node IDs (either the ancestral node of interest, which results from the
    first two above conditions; or the universal root node, resulting from the third and fourth)
    is returned as the first element of the result list.
    The remainder of the list consists of the IDs of each node in the lineage.

    This subroutine keeps passed variables to a minimum, for speed's sake; so these are globals:
    $ancestral_node; $universal_clade_node ; @parent_of ; @node_in_clade .

    XXXXAlso this could be done with a recursive subroutine call, but again I want to avoid that.XXXX No, it's recursive; otherwise it's v. messy

=cut


    my $next_node = shift;
    my @lineage = @_;

    push @lineage, $next_node;

    my $terminus = first { $_ ==  $next_node } ($ancestral_node, $universal_clade_node);

    if (defined $node_in_clade[$next_node]) { # undef means it's not known yet

        if ($node_in_clade[$next_node]) {
            # barely worthwhile sanity-check
            croak qq{conflicting clade membership for node $next_node (\$terminus == $terminus; \$node_in_clade[$next_node] == $node_in_clade[$next_node])} 
                if ((defined $terminus) && ($terminus != $ancestral_node));
            $terminus = $ancestral_node;
        }
        else {
            # barely worthwhile sanity-check
            croak qq{conflicting clade membership for node $next_node (\$terminus == $terminus; \$node_in_clade[$next_node] == $node_in_clade[$next_node])} 
                if ((defined $terminus) && ($terminus != $universal_clade_node));
            $terminus = $universal_clade_node;   
        }
    }

    if (defined $terminus) {

        my $boolean = 0;
        if ($terminus == $ancestral_node) { $boolean = 1; }

        # sets the clade membership (either known to be or known not to be) of all members of the lineage
        for my $node_id (@lineage) {

            if (defined $node_in_clade[$node_id]) {
                # sanity check
                croak qq{conflicting clade membership for node $node_id (\$terminus == $terminus; \$node_in_clade[$node_id] == $node_in_clade[$node_id])}
                    if ($node_in_clade[$node_id] != $boolean);
            }
            else {
                $node_in_clade[$node_id] = $boolean;
            }
        }

        return ($terminus, @lineage);
    }

    return get_lineage($parent_of[$next_node], @lineage);
}


my $n_iterations = 1;

print STDERR qq{\n\ncompleted first iteration; }, scalar @clade_nodes, qq{ nodes in list\n};

ITERATION:
while ($n_iterations <= $max_iterations) {

    $n_iterations++;
    my $n_added_this_iteration = 0;

    my $processed_nodes = 0;

    NODE:
    for my $node_id (0 .. $#parent_of) {

        next NODE if (!defined $parent_of[$node_id]); # some node_ids are unused in the nodes.dmp file

        my $parent_id = $parent_of[$node_id];

        if (first { $_ == $parent_id } @clade_nodes) {
            push @clade_nodes, $node_id;
            print STDERR qq{iteration $n_iterations: adding $node_id (parent is $parent_id)\n};
            $n_added_this_iteration++;
        }

        if ( !(++$processed_nodes % $report_frequency) ) {
            
            print STDERR qq{CHECKPOINT: iteration $n_iterations; $processed_nodes nodes read\n};
        }

    }

    print STDERR qq{\n\ncompleted $n_iterations iterations; },
                 qq{$n_added_this_iteration node(s) added this iteration; },
                 scalar @clade_nodes, qq{ nodes now in list\n\n\n};

    last ITERATION if (!$n_added_this_iteration);
}

print qq{# List of node $ancestral_node and all of its descendant nodes:\n},
    join(qq{\n}, @clade_nodes), qq{\n};

