"""
recon_viewer.py
View a single reconciliation using matplotlib
Credit to Justin Jiang, Trenton Wesley, and the Ran lab

Additions made by Michelle Johnson with the Bush lab

To create images, this module should be imported into interactiveAnalysis. The 
function ofamRender returns a figure object, which has the methods .save and 
.show to view the figure. The .show method will open the figure in a matplotlib 
window, and the image that gets saved is dependent on the window size (i.e., a 
call to save without a call to show may be less readable than showing the figure
and saving from the matplotlib window. Some example calls to create and save/
show a figure looks as follows:
    figure = ofamRender(genesO, originFamiliesO, speciesRtreeO, originFamilyNumber)
    figure.save("filename.extension")
or, 
    figure3099 = ofamRender(genesO, originFamiliesO, speciesRtreeO, 3099)
    figure3099.show()
    figure3099.save("ofam_plot-3099.svg")
etc.

To show a legend in the image, use
    figure = ofamRender(genesO, originFamiliesO, speciesRtreeO, originFamilyNumber, showLegend = True)

Likewise, internal labels can be turned off using
    figure = ofamRender(genesO, originFamiliesO, speciesRtreeO, originFamilyNumber, showInternalLabels = False)

There is sometimes an error where trees have cycles in them. This causes the 
object to be returned to be None.

The interface with rest of xenoGI only works for origin family objects, and is 
dependent on objects most easily accessed through interactive analysis. However,
the required objects can also be set up in an interactive workspace using the
setup function. All functions added to interface with xenoGI are located at the 
bottom of this file, under the block "NECESSARY XENOGI TRANSLATION FUNCTIONS". 

Various edits have been made to the original recon_viewer file: some fixed bugs,
others just edited the file to interface better with xenoGI inputs. Images are 
readable, but losses are no longer branches but instead nodes (red X). 
dict_to_reconciliation also returns two objects (each a separate class), one 
which holds DTS objects and the other that holds LOR objects. Most edits made
which I am uncertain of are marked with a comment containing the word "BETA". 

Current Issues:
    Legend is incorrect
    Losses currently mapped to nodes are incorrect - should be mapped to one of
        the branches under that node (which is unknown)
    X Temporal Graphs contain cycles, likely as a result of an error in build_temporal_graph
        - fixed ? by commenting out a large block

Future Improvements:
    Legend could be placed somewhere that doesn't obscure the figure
    Loss branches could actually lose said branches
    Duplications could have two branches as opposed to only one
        This would require information which is currently unavailable to us
    Text could sized appropriately (width as well as height)
"""
from __future__ import annotations
# from xenoGI.families import REARRANGEMENT

from typing import Union, Dict, Tuple, List, NamedTuple
import math

from collections import OrderedDict, namedtuple
from enum import Enum

from abc import ABC  # Abstract Base Classes

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import rcParams
from matplotlib.collections import LineCollection
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib import font_manager
import matplotlib.patheffects as PathEffects

FONT_SIZE_STRETCH = 3.0     # Parameter used in _calculate_font_size
SIGMOID_SCALE = 0.8         # Parameter used in _sigmoid

def render(host_dict: dict, parasite_dict: dict, recon_dict: dict, event_frequencies: Dict[tuple, float] = None, show_internal_labels: bool = False, 
           show_freq: bool = False, show_legend: bool = True, axes: Union[plt.Axes, None] = None):
    """ Renders a reconciliation using matplotlib
    :param host_dict:  Host tree represented in dictionary format
    :param parasite_dict:  Parasite tree represented in dictionary format
    :param recon_dict: Reconciliation represented in dictionary format
    :param event_frequencies: Dictionary that maps event tuples to their frequencies
    :param show_internal_labels: Boolean that determines wheter internal labels are shown or not
    :param show_freq: Boolean that determines whether event frequencies are shown or not
    :param axes: If specified, draw on the axes instead of creating a new figure
    :return Figure Object
    """
    host_tree, parasite_tree, consistency_type = build_trees_with_temporal_order(host_dict, parasite_dict, recon_dict)
    recon_obj, branchRecon_obj = dict_to_reconciliation(recon_dict, event_frequencies)

    # Checks to see if the trees(or reconciliation) are empty
    if host_tree is None or parasite_tree is None or recon_obj is None:
        return None

    fig = FigureWrapper(TREE_TITLE, axes)

    if show_legend:
        _create_legend(fig, consistency_type)

    # Calculates font sizes
    num_tip_nodes = len(host_tree.leaf_list()) + len(parasite_tree.leaf_list())
    font_size = _calculate_font_size(num_tip_nodes)

    root = parasite_tree.root_node
    host_lookup = host_tree.name_to_node_dict()
    parasite_lookup = parasite_tree.name_to_node_dict()

    # Populate Host Nodes with track count
    _populate_host_tracks(root, recon_obj, host_lookup)

    # Render Host Tree
    _render_host(fig, host_tree, show_internal_labels, font_size)

    # Sets the offsets between tracks on each host node
    _set_offsets(host_tree)

    # Determine the length of the longest string in the host tree's leaf list
    longest_host_name = max([len(leaf.name) for leaf in host_tree.leaf_list()])
    # Render Parasite Tree
    _render_parasite(fig, parasite_tree, recon_obj, branchRecon_obj, host_lookup, parasite_lookup, show_internal_labels, show_freq, font_size, longest_host_name)

    return fig


def _create_legend(fig: FigureWrapper, consistency_type: str):
    """
    Creates a legend on the figure
    :param fig: Figure object that visualizes trees using MatplotLib
    :param consistency_type: String that gives the consistency of a tree
    """
    fig.set_legend(LEGEND_ELEMENTS, title=consistency_type)


def _set_offsets(tree: Tree):
    """
    Populates the nodes of a Tree with an offset
    :param tree: Tree Object
    """
    pos_dict = tree.pos_dict

    for node in tree.postorder_list():
        y_1 = None
        y_0 = node.layout.row
        for logical_pos in pos_dict:
            if node.is_leaf():
                if y_0 < logical_pos[0] and node.layout.col <= logical_pos[1]:
                    if y_1 is None or y_1 > logical_pos[0]:
                        y_1 = logical_pos[0]
            else:
                y_1 = max(node.left_node.layout.row, node.right_node.layout.row)
        if y_1 is None or node.layout.node_count == 0:
            node.layout.offset = TRACK_OFFSET
        else:
            # Gives an offset based on the predicted number of horizontal tracks mapped to a host node
            # COUNT_OFFSET artificially adds extra nodes/tracks to lower the offset and pull parasite nodes closer to the host node their mapped to
            node.layout.offset = abs(y_0 - y_1) / (node.layout.node_count + COUNT_OFFSET)


def _render_host(fig: FigureWrapper, host_tree: Tree, show_internal_labels: bool, font_size: float):
    """
    Renders the host tree
    :param fig: Figure object that visualizes trees using MatplotLib
    :param host_tree: Host tree represented as a Tree object
    :param show_internal_labels: Boolean that determines whether or not the internal labels are shown
    :param font_size: Font size for the text of the tips and internal nodes of the tree
    """
    _set_host_node_layout(host_tree)
    root = host_tree.root_node
    _draw_host_handle(fig, root)
    _render_host_helper(fig, root, show_internal_labels, font_size, host_tree)


def _draw_host_handle(fig: FigureWrapper, root: Node):
    """
    Draw edge leading to root of host tree.
    :param fig: Figure object that visualizes trees using MatplotLib
    :param root: The root node of a tree object
    """
    fig.line((0, root.layout.y), (root.layout.x, root.layout.y), HOST_EDGE_COLOR)


def _render_host_helper(fig: FigureWrapper, node: Node, show_internal_labels: bool, font_size: float, host_tree: Tree):
    """
    Helper function for rendering the host tree.
    :param fig: Figure object that visualizes trees using MatplotLib
    :param node: Host Node object that will be rendered
    :param show_internal_labels: Boolean that determines whether or not the internal labels are shown
    :param font_size: Font size for the text of the tips and internal nodes of the tree
    :param host_tree: Tree object representing a Host Tree
    """
    host_tree.pos_dict[(node.layout.row, node.layout.col)] = node

    node_pos = Position(node.layout.x, node.layout.y)

    if node.is_leaf():
        text_offset = (node_pos.x + TIP_TEXT_OFFSET_X, node_pos.y)
        if node.layout.node_count == 0:
            fig.text_v2(text_offset, node.name, HOST_NODE_COLOR, size=font_size, vertical_alignment=TIP_ALIGNMENT)
        else:
            fig.text_v2(text_offset, node.name, HOST_NODE_COLOR, size=font_size, vertical_alignment=TIP_ALIGNMENT)    
    else:
        # fig.dot(node_pos, col=HOST_NODE_COLOR)  # Render host node
        fig.halfDot(node_pos, col=HOST_NODE_COLOR) # slightly smaller species node
        if show_internal_labels:
            color = transparent_color(HOST_NODE_COLOR, INTERNAL_NODE_ALPHA)
            text_xy = (node_pos.x, node_pos.y)
            fig.text_v2(text_xy, node.name, color, size=font_size*1/2, border_col=HOST_NODE_BORDER_COLOR) #the *1/2 reduces the font size
        left_x, left_y = node.left_node.layout.x, node.left_node.layout.y
        right_x, right_y = node.right_node.layout.x, node.right_node.layout.y
        fig.line(node_pos, (node_pos.x, left_y), HOST_EDGE_COLOR)
        fig.line(node_pos, (node_pos.x, right_y), HOST_EDGE_COLOR)
        fig.line((node_pos.x, left_y), (left_x, left_y), HOST_EDGE_COLOR)
        fig.line((node_pos.x, right_y), (right_x, right_y), HOST_EDGE_COLOR)
        _render_host_helper(fig, node.left_node, show_internal_labels, font_size, host_tree)
        _render_host_helper(fig, node.right_node, show_internal_labels, font_size, host_tree)


def _render_parasite(fig: FigureWrapper, parasite_tree: Tree, recon_obj: Reconciliation,  
        branch_recon: Reconciliation, host_lookup: dict, parasite_lookup: dict, 
        show_internal_labels: bool, show_freq: bool, font_size: float, longest_host_name: int):
    """
    Render the parasite tree.
    :param fig: Figure object that visualizes trees using MatplotLib
    :param parasite_tree: Parasite tree represented as a Tree object
    :param recon_obj: Reconciliation object that represents an edge-to-edge mapping from a parasite tree to a host tree
    :param branch_recon: Reconciliation object that contains 'branch' events - i.e., O R and L
    :param host_lookup: Dictionary with host node names as the key and host node objects as the values
    :param parasite_lookup: Dictionary with parasite node names as the key and parasite node objects as the values
    :param show_internal_labels: Boolean that determines whether or not the internal labels are shown
    :param show_freq: Boolean that determines wheter or not the frequencies are shown
    :param font_size: Font size for the text of the tips and internal nodes of the tree
    :param longest_host_name: The number of symbols in the longest host tree tip name
    """
    root = parasite_tree.root_node
    _render_parasite_helper(fig, root, recon_obj, branch_recon, host_lookup, parasite_lookup, show_internal_labels, show_freq, font_size, longest_host_name)


def _populate_host_tracks(node: Node, recon_obj: Reconciliation, host_lookup: dict):
    """
    :param node: Node object
    :param recon_obj: Reconciliation object that represents an edge-to-edge mapping from  a parasite tree to a host tree
    :param host_lookup: Dictionary with host node names as the key and host node objects as the values
    """
    mapping_node = recon_obj.mapping_of(node.name)
    event = recon_obj.event_of(mapping_node)
    host_name = mapping_node.host
    host_node = host_lookup[host_name]

    if not(event.event_type is EventType.DUPLICATION or event.event_type is EventType.TRANSFER):
        host_node.update_count()
    else:
        if not(_is_sharing_track(node, host_name, recon_obj)):
            host_node.update_count()

    if not(node.is_leaf()):
        _populate_host_tracks(node.left_node, recon_obj, host_lookup)
        _populate_host_tracks(node.right_node, recon_obj, host_lookup)


def _is_sharing_track(node: Node, host_name: str, recon_obj: Reconciliation):
    """
    Determines if a node is sharing it's horizontal track with its children
    :param node: Node object representing a parasite event
    :param host_name: Name of host node
    :param recon_obj: Reconciliation Object
    """
    left_host_name = recon_obj.mapping_of(node.left_node.name).host
    right_host_name = recon_obj.mapping_of(node.right_node.name).host

    return host_name == left_host_name or host_name == right_host_name


def _render_parasite_helper(fig: FigureWrapper,  node: Node, recon_obj: Reconciliation, branch_recon: Reconciliation, host_lookup: dict, parasite_lookup: dict, show_internal_labels: bool, show_freq: bool, font_size: float, longest_host_name : int):
    """
    Helper function for rendering the parasite tree.
    :param fig: Figure object that visualizes trees using MatplotLib
    :param node: Node object representing the parasite event being rendered
    :param recon_obj: Reconciliation object that represents an edge-to-edge mapping from  a parasite tree to a host tree
    :param branch_recon: Reconciliation object that only contains branch events (ORL)
    :param host_lookup: Dictionary with host node names as the key and host node objects as the values
    :param parasite_lookup: Dictionary with parasite node names as the key and parasite node objects as the values
    :param show_internal_labels: Boolean that determines whether or not the internal labels are shown
    :param show_freq: Boolean that determines wheter or not the frequencies are shown
    :param font_size: Font size for the text of the tips and internal nodes of the tree
    :param longest_host_name: The number of symbols in the longest host tree tip name
    """
    # mapping_node is of type MappingNode which associates
    # a parasite to a host in a reconciliation
    mapping_node = recon_obj.mapping_of(node.name)

    # A reconciliation has an event_of method which is an object of
    # type Event.
    event = recon_obj.event_of(mapping_node)

    host_name = mapping_node.host

    # host_lookup is a dictionary computed in the Tree class
    # that associates a host name (a string) with the correspond node
    # object for that host.  The node object contains layout information
    # which we need here.
    host_node = host_lookup[host_name]

    # Set parasite node layout
    host_row = host_node.layout.row
    # host_col = host_node.layout.col
    # host_x = host_node.layout.x
    host_y = host_node.layout.y
    node.set_layout(row=host_row, x=node.layout.col, y=host_y)

    # Render parasite node and recurse if not a leaf
    if node.is_leaf():
        node.layout.y += host_node.get_and_update_track(Track.HORIZONTAL) * host_node.layout.offset
        _render_parasite_node(fig, node, event, font_size, longest_host_name)
        # draw branch events - beta
        if node.name in branch_recon._parasite_map:
            mapping_nodes = branch_recon.mapping_of(node.name)
            for mapping_node in mapping_nodes: # MICHELLE (changed _parasite_map to form (str, list))
                branch_event = branch_recon.event_of(mapping_node)
                _render_branch_event(node, mapping_node, branch_event, host_lookup, fig)
        return

    # If the Node is in their own track, change their position
    if not(_is_sharing_track(node, host_name, recon_obj)):
        node.layout.y += host_node.layout.h_track * host_node.layout.offset

    left_node, right_node = _get_children(node, recon_obj, parasite_lookup)

    _render_parasite_helper(fig, left_node, recon_obj, branch_recon, host_lookup,
        parasite_lookup, show_internal_labels, show_freq, font_size, longest_host_name)
    _render_parasite_helper(fig, right_node, recon_obj, branch_recon, host_lookup,
        parasite_lookup, show_internal_labels, show_freq, font_size, longest_host_name)

    # Checking to see if left node is mapped to the same host node as parent
    if node.layout.row == left_node.layout.row:
        node.set_layout(y=left_node.layout.y)
    elif node.layout.row == right_node.layout.row:
        node.set_layout(y=right_node.layout.y)
    elif event.event_type is EventType.TRANSFER:
        node.layout.y = host_node.layout.y + host_node.layout.h_track * host_node.layout.offset

    #Checks to see if transfer node is inconsistent and if it can be fixed
    if event.event_type is EventType.TRANSFER:
        min_col = host_lookup[recon_obj.mapping_of(right_node.name).host].parent_node.layout.col
        if min_col >= node.layout.col:
            _fix_transfer(node, left_node, right_node, host_node, host_lookup, parasite_lookup, recon_obj)

    _render_parasite_branches(fig, node, recon_obj, branch_recon, host_lookup, parasite_lookup)
    _render_parasite_node(fig, node, event, font_size, longest_host_name, show_internal_labels, show_freq)


def _fix_transfer(node: Node, left_node, right_node: Node, host_node: Node, host_lookup: dict, parasite_lookup: dict, recon_obj: Reconciliation, node_col: float = None, offset_number: int = 1):
    """
    Checks to see in tranfer node is inconsistent and the tries to fix node if it can be slid down the host edge
    The tries to push a given node forward if possible to correct the assumed inconsistency
    :param node: Node object representing the parasite event being rendered
    :param node: Right node of the node object
    :param left_node: Left node of the node object
    :param right_node: Right node of the node object
    :param host_node: Node object represeting a host that the parasite node is mapped to
    :param host_lookup: Dictionary with host node names as the key and host node objects as the values
    :param parasite_lookup: Dictionary with parasite node names as the key and parasite node objects as the values
    :param recon_obj: Reconciliation object that represents an edge-to-edge mapping from  a parasite tree to a host tree
    :param node_col: Column of the node, used when function is called recursively to check the next available transfer spot
    :param offset_number: Used to push node to an available transfer spot
    """
    min_col = host_lookup[recon_obj.mapping_of(right_node.name).host].parent_node.layout.col
    max_col = host_node.layout.col
    node_col = node.layout.col
    max_col = min(host_node.layout.col, left_node.layout.col)
    if not(node_col):
        node_col = node.layout.col

    # Checks to see if transfer is inconsistent and if the inconsistency can be fixed by sliding the transfer node down the host edge
    if min_col >= node_col and min_col < max_col and not(_is_sharing_track(node, host_node.name, recon_obj)):
        node.set_layout(col=min_col+0.5, x=min_col+0.5)
    if min_col < max_col:
        new_value = min_col + PUSHED_NODE_OFFSET * offset_number
        if _is_col_taken(new_value, host_lookup, parasite_lookup):
            _fix_transfer(node, left_node, right_node, host_node, host_lookup, parasite_lookup, recon_obj, node_col=new_value, offset_number = offset_number + 1)
        else:
            node.set_layout(col=new_value, x=new_value)


def _is_col_taken(node_col, host_lookup, parasite_lookup):
    """
    Checks to see if a node is already in a given col
    :param node_col: Column of a given node
    :param host_lookup: Dictionary with host node names as the key and host node objects as the values
    :param parasite_lookup: Dictionary with parasite node names as the key and parasite node objects as the values
    :return True if there is already a node at the given column and False otherwise
    """
    for key in host_lookup:
        if host_lookup[key].layout.col == node_col:
            return True
    for key in parasite_lookup:
        if parasite_lookup[key].layout.col == node_col:
            return True
    return False


def _render_parasite_node(fig: FigureWrapper,  node: Node, event: Event, font_size: float, longest_host_name : int, show_internal_labels: bool = False, show_freq: bool = False):
    """
    Renders a single parasite node
    :param fig: Figure object that visualizes trees using MatplotLib
    :param node: parasite Node that will be rendered
    :param event: Event object that gives event for the given Node
    :param font_size: Font size for text
    :param show_internal_labels: Boolean that determines whether or not the internal labels are shown
    :param show_freq: Boolean that determines wheter or not the frequencies are shown
    :param longest_host_name: The number of symbols in the longest host tree tip name
    """

    node_pos = Position(node.layout.x, node.layout.y)
    render_color, render_shape = _event_color_shape(event)


    if node.is_leaf():
        fig.text_v2((node.layout.x + TIP_TEXT_OFFSET_X, node.layout.y), "-"*(3+longest_host_name)+node.name, render_color, size=font_size, vertical_alignment=TIP_ALIGNMENT)
        return

    fig.dot(node_pos, col=render_color, marker=render_shape)
    text = ''
    text_color = transparent_color(render_color, INTERNAL_NODE_ALPHA)
    if show_internal_labels and show_freq:
        text = node.name + ', ' + _get_frequency_text(event.freq)
    elif show_internal_labels:
        text = node.name
    elif show_freq:
        if event.freq:
            text = _get_frequency_text(event.freq)
        else:
            raise RuntimeError("Could not render reconciliation: show_freq is True but event.freq is None")
    if text:
        fig.text_v2(node_pos, text, text_color, size=font_size * (1/2), border_col=PARASITE_NODE_BORDER_COLOR) # *1/2 reduces the font size


def _get_frequency_text(frequency: float):
    """
    Give the frequency as a string in percentage form
    :param frequency: The frequency of an event
    :return a string that has the frequency as a percentage
    """
    return str(round(frequency * 100))


def _calculate_font_size(num_tip_nodes: int):
    """
    Calculates the font_size
    :param num_tip_nodes: Number of tip nodes in a tree
    :return the font size for the tips and internal nodes of a tree
    """
    x = (START_SIZE - num_tip_nodes) / STEP_SIZE
    return FONT_SIZE_STRETCH * _sigmoid(x)  # 3.0 is a magic value that can be adjusted


def _sigmoid(x: float):
    """
    sigmoid Function
    :param x: A number to be plugged into the function
    :return a sigmoid value based on the input value, x
    """
    return (1 / (1 + math.e**(-SIGMOID_SCALE*x)))  # 0.8 is a magic value that can be adjusted


def _render_parasite_branches(fig: FigureWrapper,  node: Node, recon_obj: Reconciliation, branch_recon: Reconciliation, host_lookup: dict, parasite_lookup: dict):
    """
    Very basic branch drawing
    :param fig: Figure object that visualizes trees using MatplotLib
    :param node: Node object representing the parasite event being rendered
    :param recon_obj: Reconciliation object that represents an edge-to-edge mapping from  a parasite tree to a host tree
    :param branch_recon: Reconciliation object representing ORL events (branch events)
    :param host_lookup: Dictionary with host node names as the key and host node objects as the values
    :param parasite_lookup: Dictionary with parasite node names as the key and parasite node objects as the values
    """
    node_pos = Position(node.layout.x, node.layout.y)

    left_node, right_node = _get_children(node, recon_obj, parasite_lookup)

    right_pos = Position(right_node.layout.x, right_node.layout.y)

    mapping_node = recon_obj.mapping_of(node.name)
    event = recon_obj.event_of(mapping_node)

    if event.event_type is EventType.COSPECIATION:
        _render_cospeciation_branch(node, host_lookup, parasite_lookup, recon_obj, fig)
    elif event.event_type is EventType.DUPLICATION:
        _render_cospeciation_branch(node, host_lookup, parasite_lookup, recon_obj, fig) # changed to have two children for duplication, originally # _connect_children(node, host_lookup, parasite_lookup, recon_obj, fig)
    elif event.event_type is EventType.TRANSFER:
        left_pos = Position(left_node.layout.x, left_node.layout.y)
        # determine which side is the same, draw the other side as transfer
        if node_pos.y != right_pos.y:
            _render_transfer_branch(node_pos, right_pos, fig, node, host_lookup, recon_obj, right_node)
            _render_normal_branch(node_pos, left_pos, fig) # fix issues where the transfer branch tracks all the way back to origin instead of desired stop
        else:
            _render_transfer_branch(node_pos, left_pos, fig, node, host_lookup, recon_obj, right_node)
            _render_normal_branch(node_pos, right_pos, fig)
        # _connect_child_to_parent(node, left_node, host_lookup, recon_obj, fig) # why do I need this?? MICHELLE
    else:
        raise ValueError('%s is not an known event type' % event.event_type)
    # Losses are implicitly implied and are not mapped to any specific node

    # draw branch events - beta
    if node.name in branch_recon._parasite_map:
        mapping_nodes = branch_recon.mapping_of(node.name)
        for mapping_node in mapping_nodes:
            branch_event = branch_recon.event_of(mapping_node)
            _render_branch_event(node, mapping_node, branch_event, host_lookup, fig)

def _render_branch_event(parasite_node, mapping_node, event, host_lookup, fig):
    """
    Draws a node on the branch leading to the given parasite_node.
    :param parasite_node:  the node the event happens on/before
    :param mapping_node:  node mapping the parasite_node to a species_node
    :param event:  the event type
    :param host_lookup:  dictionary with host node names as key and host node objects as values
    :param fig:  Figure Object
    """
    # node is a PARASITE node
    host_name = mapping_node.host
    host_node = host_lookup[host_name]

    parent_node = None

    if event.event_type == EventType.ORIGIN:
        parent_node = host_node.parent_node
        # ERROR CAN OCCUR IF HOST_NODE IS s0
        if parent_node == None:
            offset = host_node.layout.x - parasite_node.layout.x
            parent_node = Node('parentOfOrigin')
            parent_node.layout = NodeLayout()
            parent_node.set_layout(x=parasite_node.layout.x - (offset), \
                                    y=parasite_node.layout.y)
            
    else:
        parent_node = parasite_node.parent_node
    
    host_parent_pos = Position(parent_node.layout.x, \
            parent_node.layout.y)

    '''try:
        # parent_node = parasite_node.parent_node # fixing position ??
        host_parent_pos = Position(parent_node.layout.x, \
            parent_node.layout.y)
        print('try: ' + str(parasite_node))
    except: # origin has no parent ?
        offset = host_node.layout.offset
        host_parent_pos = Position(host_node.layout.x - offset, \
            host_node.layout.y)
        print('except: ' + str(parasite_node))'''

    current_pos = Position(parasite_node.layout.x, parasite_node.layout.y)

    # node_pos = Position(host_node.layout.x, host_node.layout.y)
    # mid_pos = Position(node_pos.x, current_pos.y)

    # get rearrangement nodes, do not have loss nodes
    if event.event_type is not EventType.LOSS:
        render_color, render_shape = _event_color_shape(event)
        mid_pos = Position(host_parent_pos.x, current_pos.y)
        fig.dot(mid_pos, col=render_color, marker=render_shape)

    if event.event_type is EventType.ORIGIN:
        fig.line(mid_pos, current_pos, PARASITE_EDGE_COLOR)

def _connect_children(node: Node, host_lookup: dict, parasite_lookup: dict, recon_obj: Reconciliation, fig: FigureWrapper):
    """
    Connects the children of a node
    :param node: Node object representing a parasite event
    :param host_lookup: Dictionary with host node names as the key and host node objects as the values
    :param recon_obj: Reconciliation object that represents an edge-to-edge mapping from  a parasite tree to a host tree
    :param fig: Figure object that visualizes trees using MatplotLib
    """
    print(node)
    left_node, right_node = _get_children(node, recon_obj, parasite_lookup)
    print(left_node)
    print(right_node)
    # FUNKY EDITING, JUST TO SEE WHAT HAPPENS
    print(node.layout.y)
    node.layout.y += TRACK_OFFSET / 2
    print(node.layout.y)
    _connect_child_to_parent(node, left_node, host_lookup, recon_obj, fig)
    node.layout.y += (TRACK_OFFSET * -1)
    print(node.layout.y)
    _connect_child_to_parent(node, right_node, host_lookup, recon_obj, fig)


def _render_loss_branch(node_pos: Position, next_pos: Position, fig: FigureWrapper):
    """
    Renders a loss branch given a two positions
    :param node_pos: x and y position of a node
    :param next_pos: x and y position of another node
    :param fig: Figure object that visualizes trees using MatplotLib
    """
    # Create vertical line to next node
    mid_pos = Position(node_pos.x, next_pos.y)
    fig.line(node_pos, mid_pos, LOSS_EDGE_COLOR, linestyle='--')
    fig.line(mid_pos, next_pos, PARASITE_EDGE_COLOR)

def _render_normal_branch(node_pos: Position, next_pos:Position, fig: FigureWrapper):
    # Create vertical line to next node
    mid_pos = Position(node_pos.x, next_pos.y)
    fig.line(node_pos, mid_pos, PARASITE_EDGE_COLOR) 
    fig.line(mid_pos, next_pos, PARASITE_EDGE_COLOR)

def _render_cospeciation_branch(node: Node, host_lookup: dict, parasite_lookup: dict, recon_obj: Reconciliation, fig: FigureWrapper):
    """
    Renders the cospeciation branch.
    :param node: Node object representing the parasite event being rendered
    :param host_lookup: Dictionary with host node names as the key and host node objects as the values
    :param parasite_lookup: Dictionary with parasite node names as the key and parasite node objects as the values
    :param recon_obj: Reconciliation object that represents an edge-to-edge mapping from  a parasite tree to a host tree
    :param fig: Figure object that visualizes trees using MatplotLib
    """
    left_node, right_node = _get_children(node, recon_obj, parasite_lookup)

    node_pos = Position(node.layout.x, node.layout.y)
    left_pos = Position(left_node.layout.x, left_node.layout.y)
    right_pos = Position(right_node.layout.x, right_node.layout.y)

    mapping_node = recon_obj.mapping_of(node.name)
    host_node = host_lookup[mapping_node.host]

    # Update h_track
    host_node.get_and_update_track(Track.HORIZONTAL)

    left_mapping_node = recon_obj.mapping_of(left_node.name)
    left_host_node = host_lookup[left_mapping_node.host]

    right_mapping_node = recon_obj.mapping_of(right_node.name)
    right_host_node = host_lookup[right_mapping_node.host]
    # Draw left node
    offset = host_node.layout.offset
    if host_node.left_node.name == left_host_node.name:
        _render_curved_line_to(node_pos, left_pos, fig)
        if host_node.layout.lower_v_track < (host_node.layout.x - node_pos.x) / offset:
            host_node.layout.lower_v_track += (host_node.layout.x - node_pos.x) / offset + offset
    else:
        # stop_row = host_node.left_node.layout.row # -beta test, parse in left & right as stop nodes to catch cases where left/right labels are swapped
        stop_rows = (host_node.left_node.layout.row, host_node.right_node.layout.row)
        _connect_child_to_parent(node, left_node, host_lookup, recon_obj, fig, stop_rows=stop_rows) # -beta test

    # Draw Right node
    if host_node.right_node.name == right_host_node.name:
        _render_curved_line_to(node_pos, right_pos, fig)
        if host_node.layout.upper_v_track < (host_node.layout.x - node_pos.x) / offset:
            host_node.layout.upper_v_track += (host_node.layout.x - node_pos.x) / offset + offset
    else:
        # stop_row = host_node.right_node.layout.row # -beta test, same as above
        stop_rows = (host_node.left_node.layout.row, host_node.right_node.layout.row)
        _connect_child_to_parent(node, right_node, host_lookup, recon_obj, fig, stop_rows=stop_rows) # originally stop_row


def _get_children(node: Node, recon_obj: Reconciliation, parasite_lookup: dict):
    """
    Gets the children of a node in the order they appear in the mapping node.
    :param node: Node object representing a parasite event
    :param recon_obj: Reconciliation Object
    :param parasite_lookup: Dictionary with parasite node names as the key and parasite node objects as the values
    :return A tuple consisting of the left node and right node
    """
    mapping_node = recon_obj.mapping_of(node.name)
    event = recon_obj.event_of(mapping_node)
    left_mapping_node = event.left
    right_mapping_node = event.right
    left_node_name = left_mapping_node.parasite
    right_node_name = right_mapping_node.parasite

    left_node = parasite_lookup[left_node_name]
    right_node = parasite_lookup[right_node_name]

    return left_node, right_node


def _render_curved_line_to(node_pos: Position, other_pos: Position, fig: FigureWrapper):
    """
    Renders a curved line from one point to another
    :param node_pos: x and y position of a node
    :param other_pos: x and y position of another node
    :param fig: Figure object that visualizes trees using MatplotLib
    """
    mid_pos = Position(node_pos.x, other_pos.y)
    fig.line(node_pos, mid_pos, PARASITE_EDGE_COLOR)
    fig.line(mid_pos, other_pos, PARASITE_EDGE_COLOR)


def _render_transfer_branch(node_pos: Position, right_pos: Position, fig: FigureWrapper,  node: Node, host_lookup: dict, recon_obj: Reconciliation, right_node: Node):
    """
    Renders a transfer branch
    :param node_xy: x and y position of a node
    :param right_pos: x and y position of the right child of a node
    :param fig: Figure object that visualizes trees using MatplotLib
    :param node: Node object representing the parasite event being rendered
    :param host_lookup: Dictionary with host node names as the key and host node objects as the values
    :param recon_obj: Reconciliation object that represents an edge-to-edge mapping from  a parasite tree to a host tree
    :param right_node: The right node object of node
    """

    child_mapping_node = recon_obj.mapping_of(right_node.name)
    child_host_node = host_lookup[child_mapping_node.host]

    # Check temporal consistency of transfer event
    if child_host_node.parent_node.layout.col < node.layout.col:
        # Draw right node, which is transfered
        mid_pos = Position(node_pos.x, right_pos.y)  # xy coords of midpoint
        y_midpoint = abs(mid_pos.y + node_pos.y) / 2  # value of midpoint between mid_xy and parent node

        # Determine if transfer is upwards or downwards, and draw triangle accordingly
        is_upwards = True if y_midpoint < mid_pos.y else False
        arrow_pos = Position(node_pos.x, y_midpoint)
        if is_upwards:
            fig.triangle(arrow_pos, TRANSFER_NODE_COLOR)
        else:
            fig.triangle(arrow_pos, TRANSFER_NODE_COLOR, rotation=DOWN_ARROW_ROTATION)

        # Draw branch to midpoint, then draw branch to child
        fig.line(node_pos, mid_pos, TRANSFER_NODE_COLOR)
        fig.line(mid_pos, right_pos, TRANSFER_NODE_COLOR)
    else:
        transfer_edge_color = MAROON # previously: transparent_color(PARASITE_EDGE_COLOR, TRANSFER_TRANSPARENCY) -beta
        # note: this appears to have made no difference on plots ??
        fig.arrow_segment(node_pos, right_pos, transfer_edge_color)
        fig.line(node_pos, right_pos, transfer_edge_color)


def _connect_child_to_parent(node: Node, child_node: Node, host_lookup: dict, recon_obj: Reconciliation, fig: FigureWrapper, stop_rows: list = []): # stop_row: float = None): MICHELLE
    """
    Connects a child node to its parent node
    :param node: Node object representing a parasite event
    :param child_node: The child node object of a given node
    :param host_lookup: Dictionary with host node names as the key and host node objects as the values
    :param recon_obj: Reconciliation object that represents an edge-to-edge mapping from  a parasite tree to a host tree
    :param fig: Figure object that visualizes trees using MatplotLib
    :param stop_row: row number to stop line drawing on
    """
    # get map between node (parasite) and host
    mapping_node = recon_obj.mapping_of(child_node.name)
    host_node = host_lookup[mapping_node.host]
      
    # switch stop rows from singular to multiple, to catch left/right issues in cospeciation events
    # if no stop_row specified, use input node as stopping point
    #if stop_row == None:
    if stop_rows == []:
        stop_rows = [node.layout.row] # originally singular
    
    # start at the child
    current_pos = Position(child_node.layout.x, child_node.layout.y)

    # """ REMOVAL START
    # if the host we are looking at is not in the same row as our given node, we loop
    # while host_node.layout.row != stop_row and host_node.parent_node:
    while host_node.layout.row not in stop_rows and host_node.parent_node:
        parent_node = host_node.parent_node
        # if the parent node is below the current node, we're on the upper track
        if parent_node.layout.row < host_node.layout.row:
            v_track = parent_node.get_and_update_track(Track.UPPER_VERTICAL)
        # else, we are on the lower track
        else:
            v_track = parent_node.get_and_update_track(Track.LOWER_VERTICAL)
            # update track until we are where the upper track is (???)
            while v_track < parent_node.layout.upper_v_track:
                v_track = parent_node.get_and_update_track(Track.LOWER_VERTICAL)
        # get the horizontal track
        h_track = parent_node.get_and_update_track(Track.HORIZONTAL)
        # create offset!
        offset = parent_node.layout.offset

        # get the position to draw to, should be next to the species node before the current node
        # sub_parent_pos = Position(parent_node.layout.x - (offset * v_track), \
        #     parent_node.layout.y + (offset * h_track))
        
        sub_parent_pos = Position(parent_node.layout.x - offset, parent_node.layout.y + offset) # commented out, written by Michelle as a test (beta)

        _render_loss_branch(sub_parent_pos, current_pos, fig) # should be _loss_branch MICHELLE

        # iterate such that we start at the last host node we drew to
        host_node = parent_node
        current_pos = sub_parent_pos # """ REMOVAL END
    
    node_pos = Position(node.layout.x, node.layout.y)
    mid_pos = Position(node_pos.x, current_pos.y)

    fig.line(node_pos, mid_pos, PARASITE_EDGE_COLOR)
    fig.line(mid_pos, current_pos, PARASITE_EDGE_COLOR)


def _event_color_shape(event: Event):
    """
    Gives the color and shape for drawing event, depending on event type
    :param event: Event object
    :return A tuple with the color and shape of an event
    """
    if event.event_type is EventType.TIPTIP:
        return LEAF_NODE_COLOR, LEAF_NODE_SHAPE
    if event.event_type is EventType.COSPECIATION:
        return COSPECIATION_NODE_COLOR, COSPECIATION_NODE_SHAPE
    if event.event_type is EventType.DUPLICATION:
        return DUPLICATION_NODE_COLOR, DUPLICATION_NODE_SHAPE
    if event.event_type is EventType.TRANSFER:
        return TRANSFER_NODE_COLOR, TRANSFER_NODE_SHAPE
    if event.event_type is EventType.ORIGIN:
        return ORIGIN_NODE_COLOR, ORIGIN_NODE_SHAPE
    if event.event_type is EventType.REARRANGEMENT:
        return REARRANGEMENT_NODE_COLOR, REARRANGEMENT_NODE_SHAPE
    if event.event_type is EventType.LOSS:
        return LOSS_COLOR_SHAPE
    return None, None 


def _set_host_node_layout(host_tree: Tree):
    """
    Sets the logicalRow and logicalCol values of each Node in host_tree.
    Assumes that each host Node has its order set already and this function
    uses those values and structure of the tree to set the logicalRow and logicalCol
    :param host_tree:  A Tree object representing the host tree
    :return None
    """
    # Sets logical row values for leaves in the order they appear in the list of host tree leaves
    logical_row_counter = 0
    for leaf in host_tree.leaf_list():
        leaf.layout.row = logical_row_counter
        leaf.layout.x = leaf.layout.col  # This can be scaled if desired
        leaf.layout.y = leaf.layout.row  # This can be scaled if desired
        logical_row_counter += 1
    # Helper function to assign row values, postorder traversal
    _set_internal_host_nodes(host_tree.root_node)


def _set_internal_host_nodes(node: Node):
    """
    Helper function for set_host_node_layout
    :param node: Host Node object that will be rendered
    """
    if node.is_leaf():
        return
    _set_internal_host_nodes(node.left_node)
    _set_internal_host_nodes(node.right_node)
    node.layout.row = (node.left_node.layout.row + node.right_node.layout.row) / 2
    node.layout.x = node.layout.col  # This can be scaled if desired
    node.layout.y = node.layout.row  # This can be scaled if desired

##################################################################
##################################################################
#################$ NECESSARY UTILS FUNCTIONS #####################
##################################################################
##################################################################

#__all__ = ['dict_to_tree', 'dict_to_reconciliation', 'build_trees_with_temporal_order']

class ConsistencyType(Enum):
    """ Defines type of the temporal consistency of a reconciliation """
    STRONG_CONSISTENCY = 1
    WEAK_CONSISTENCY = 2
    NO_CONSISTENCY = 3

    def __str__(self):
        if self == ConsistencyType.STRONG_CONSISTENCY:
            return "Strong time-consistent"
        elif self == ConsistencyType.WEAK_CONSISTENCY:
            return "Weak time-consistent"
        else:
            return "Not time-consistent"

# Utility functions that covert from dictionaries to objects


def dict_to_tree(tree_dict: dict, tree_type: TreeType) -> Tree:
    """
    :param tree_dict: An edge-based representation of a tree as in the example above.
    :param tree_type: tree.TreeType.{HOST, PARASITE} indicating the type of the tree.
        This is used to determine if the handle of the tree is "hTop" (host) or
        "pTop" (parasite)
    :return: A representation of the tree in Tree format (see tree.py)
    """
    # pull the root from the tree dictionary
    root = ''
    for geneKey in tree_dict:
        parent_name = tree_dict[geneKey][0]
        # root will have a parent named 'h_root' or 'p_root'
        if parent_name == 'h_root' or parent_name == 'p_root':
            root = geneKey
    output_tree = Tree()
    output_tree.tree_type = tree_type
    output_tree.root_node = _dict_to_tree_helper(tree_dict, root)
    return output_tree


def _dict_to_tree_helper(tree_dict, root_edge):
    """
    Helper function for dict_to_tree.
    """
    root_name = tree_dict[root_edge][1]
    new_node = Node(root_name)
    left_edge = tree_dict[root_edge][2]
    right_edge = tree_dict[root_edge][3]

    if left_edge is None and right_edge is None:
        return new_node
    else:
        new_left_node = _dict_to_tree_helper(tree_dict, left_edge)
        new_right_node = _dict_to_tree_helper(tree_dict, right_edge)
        new_node.left_node = new_left_node
        new_left_node.parent_node = new_node
        new_node.right_node = new_right_node
        new_right_node.parent_node = new_node
        return new_node

# ReconGraph utilities


def _find_roots(old_recon_graph) -> List[MappingNode]:
    not_roots = set()
    for mapping in old_recon_graph:
        for event in old_recon_graph[mapping]:
            etype, left, right = event
            if etype in 'SDTO':
                not_roots.add(left)
                not_roots.add(right)
            elif etype == 'L' or etype == 'R':
                child = left
                not_roots.add(child)
            elif etype == 'C':
                pass
            else:
                raise ValueError('%s not in "SDTLC' % etype)
    roots = []
    for mapping in old_recon_graph:
        if mapping not in not_roots:
            roots.append(mapping)
    return roots


def dict_to_reconciliation(old_recon: Dict[Tuple, List], event_frequencies: Dict[tuple, float] = None):
    """
    Convert the old reconciliation graph format to Reconciliation.
    Example of old format:
    old_recon = {
        ('n0', 'm2'): [('S', ('n2', 'm3'), ('n1', 'm4'))],
        ('n1', 'm4'): [('C', (None, None), (None, None))],
        ('n2', 'm3'): [('T', ('n3', 'm3'), ('n4', 'm1'))],
        ('n3', 'm3'): [('C', (None, None), (None, None))],
        ('n4', 'm1'): [('C', (None, None), (None, None))],
    }
    Returns recon_obj (SDTC) and branchRecon_obj (ORL)
    """
    # check that we have only one root
    roots = _find_roots(old_recon)
    if len(roots) > 1:
        raise ValueError("old_recon has many roots")
    root = roots[0]
    # initialize a normal reconciliation object and one for branch events
    recon_obj = Reconciliation(root)
    branchRecon_obj = BranchReconciliation(root)
    event = None
    for mapping in old_recon:
        host, parasite = mapping
        # iterate through all events
        for event_tuple in old_recon[mapping]:
            etype, left, right = event_tuple
            mapping_node = MappingNode(host, parasite)
            # events with two children
            if etype in 'SDTO':
                left_parasite, left_host = left
                right_parasite, right_host = right
                left_mapping = MappingNode(left_parasite, left_host)
                right_mapping = MappingNode(right_parasite, right_host)
                if etype == 'S':
                    event = Cospeciation(left_mapping, right_mapping)
                if etype == 'D':
                    event = Duplication(left_mapping, right_mapping)
                if etype == 'T':
                    event = Transfer(left_mapping, right_mapping)
                if etype == 'O':
                    event = Origin(left_mapping, right_mapping)
            # events with one child
            elif etype in 'LR':
                child_parasite, child_host = left
                child_mapping = MappingNode(child_parasite, child_host)
                if etype == 'L':
                    event = Loss(child_mapping)
                if etype == 'R':
                    event = Rearrangement(child_mapping)
            # leaf events
            elif etype == 'C':
                event = TipTip()
            else:
                raise ValueError('%s not in "SDTLCOR"' % etype)
            if event_frequencies is not None:
                event._freq = event_frequencies[event_tuple]
            # add event to matching reconciliation object
            if etype in 'SDTC':
                recon_obj.set_event(mapping_node, event)
            elif etype in 'ORL':
                branchRecon_obj.set_event(mapping_node, event)
    return recon_obj, branchRecon_obj

# Temporal ordering utilities


def build_trees_with_temporal_order(host_dict: dict, parasite_dict: dict, reconciliation: dict) \
        -> Tuple[Tree, Tree, ConsistencyType]:
    """
    This function uses topological sort to order the nodes inside host and parasite tree.
    The output trees can be used for visualization.

    :param host_dict: host tree dictionary
    :param parasite_dict: parasite tree dictionary
    :param reconciliation: reconciliation dictionary
    :return: a Tree object of type HOST, a Tree object of type PARASITE, and the
             ConsistencyType of the reconciliation. If the reconciliation is either
             strongly or weakly temporally consistent, then the tree objects have their
             nodes populated and the nodes will contain the temporal order information
             in the layout field. If the reconciliation is not temporally consistent, the
             function returns None, None, ConsistencyType.NO_CONSISTENCY.
    """

    consistency_type = ConsistencyType.NO_CONSISTENCY

    # find the temporal order for host and parasite nodes using strong temporal constraints
    temporal_graph = build_temporal_graph(host_dict, parasite_dict, reconciliation)
    ordering_dict = topological_order(temporal_graph)

    if ordering_dict is not None:
        consistency_type = ConsistencyType.STRONG_CONSISTENCY
    else:
        # the reconciliation is not strongly consistent, we relax the temporal constraints
        temporal_graph = build_temporal_graph(host_dict, parasite_dict, reconciliation, False)
        ordering_dict = topological_order(temporal_graph)
        if ordering_dict is not None:
            consistency_type = ConsistencyType.WEAK_CONSISTENCY

    host_tree_object = dict_to_tree(host_dict, TreeType.HOST)
    parasite_tree_object = dict_to_tree(parasite_dict, TreeType.PARASITE)

    # if there is a valid temporal ordering, we populate the layout with the order corresponding to the node
    if consistency_type != ConsistencyType.NO_CONSISTENCY:
        # calculate the temporal order for leaves, which all have the largest order
        max_order = 1
        for node in ordering_dict:
            if max_order < ordering_dict[node]:
                max_order = ordering_dict[node]
        leaf_order = max_order + 1
        populate_nodes_with_order(host_tree_object.root_node, TreeType.HOST, ordering_dict, leaf_order)
        populate_nodes_with_order(parasite_tree_object.root_node, TreeType.PARASITE, ordering_dict, leaf_order)
        return host_tree_object, parasite_tree_object, consistency_type
    else:
        return None, None, ConsistencyType.NO_CONSISTENCY


def create_parent_dict(host_dict: dict, parasite_dict: dict):
    """
    :param host_dict:  host tree dictionary
    :param parasite_dict:  parasite tree dictionary
    :return: A dictionary that maps the name of a child node to the name of its parent
             for both the host tree and the parasite tree.
    """
    parent_dict = {}
    for edge_name in host_dict:
        child_node = _bottom_node(host_dict[edge_name])
        parent_node = _top_node(host_dict[edge_name])
        parent_dict[child_node] = parent_node
    for edge_name in parasite_dict:
        child_node = _bottom_node(parasite_dict[edge_name])
        parent_node = _top_node(parasite_dict[edge_name])
        parent_dict[child_node] = parent_node
    return parent_dict


def build_formatted_tree(tree_dict):
    """
    :param tree_dict:  a tree dictionary
    :return: A temporal graph that contains all the temporal relations implied by
             the tree. Each key is a node tuple of the form (name, type) where name
             is a string representing the name of a parasite or host tree INTERNAL 
             node and type is either TreeType.HOST or TreeType.PARASITE which are 
             defined in py. The associated value is a list of node tuples that
             are the children of this node tuple in the tree.
    """
    tree_type = None
    # look for 'p_root' as a parent. If found, set treeType to parasite
    for key in tree_dict:
        if tree_dict[key][0] == 'p_root':
            tree_type = TreeType.PARASITE
    if tree_type == None:
        tree_type = TreeType.HOST

    formatted_tree = {}
    for edge_name in tree_dict:
        edge_four_tuple = tree_dict[edge_name]
        # the temporal graph does not contain leaves as keys
        if _is_leaf_edge(edge_four_tuple):
            continue
        # the temporal graph contains internal node tuples as keys,
        # and their children nodes tuples as values
        node_name = _bottom_node(edge_four_tuple)
        left_child_name = edge_four_tuple[2]
        right_child_name = edge_four_tuple[3]
        formatted_tree[(node_name, tree_type)] = [(left_child_name, tree_type), (right_child_name, tree_type)]
    return formatted_tree


def uniquify(elements):
    """
    :param elements:  a list whose elements might not be unique
    :return: A list that contains only the unique elements of the input list. 
    """
    return list(set(elements))


def build_temporal_graph(host_dict: dict, parasite_dict: dict, reconciliation: dict, add_strong_constraints = True):
    """
    :param host_dict:  host tree dictionary
    :param parasite_dict:  parasite tree dictionary
    :param reconciliation:  reconciliation dictionary
    :param add_strong_constraints:  a boolean indicating whether we are using the strongest
                                    temporal constraints, i.e. adding the constraints implied
                                    by a transfer event
    :return: The temporal graph which is defined as follows:
        Each key is a node tuple of the form (name, type) where name is a string representing
        the name of a parasite or host tree INTERNAL node and type is either TreeType.HOST or 
        TreeType.PARASITE which are defined in py. 
        Note that leaves of the host and parasite trees are not considered here.
        The associated value is a list of node tuples that are the children of this node tuple
        in the temporal graph.
    
    There's currently a huge can of worms in this function, certain nodes are being added when
    they shouldn't be, hence all the nasty print statements.
    """
    # create a dictionary that maps each host and parasite node to its parent
    parent = create_parent_dict(host_dict, parasite_dict)
    # create temporal graphs for the host and parasite tree
    temporal_host_tree = build_formatted_tree(host_dict)
    temporal_parasite_tree = build_formatted_tree(parasite_dict)
    # initialize the final temporal graph to the combined temporal graphs of host and parasite tree
    temporal_graph = temporal_host_tree
    temporal_graph.update(temporal_parasite_tree)
    # add temporal relations implied by each node mapping and the corresponding event
    for node_mapping in reconciliation:
        parasite, host = node_mapping
        host_parent = parent[host]
        # get the event corresponding to this node mapping
        event_tuple = reconciliation[node_mapping][0]
        event_type = event_tuple[0]
        # loss-type events are included under new two-reconciliations model
        # if the node_mapping is not a leaf_mapping, we add the first relation
        # if event_type != 'C':
        if event_type not in 'CL':
            if (parasite, TreeType.PARASITE) in temporal_graph:
                temporal_graph[(parasite, TreeType.PARASITE)].append((host, TreeType.HOST))
            else:
                temporal_graph[(parasite, TreeType.PARASITE)] = [(host, TreeType.HOST)]
            '''if (host, TreeType.HOST) in temporal_graph:
                if (parasite, TreeType.PARASITE) not in temporal_graph[(host, TreeType.HOST)]:
                    if (parasite, TreeType.PARASITE) in temporal_graph:
                        print("1")
                        temporal_graph[(parasite, TreeType.PARASITE)].append((host, TreeType.HOST))
                    else: # add leaves specifically for RL events on branches that lead to leaves
                        print("2")
                        temporal_graph[(parasite, TreeType.PARASITE)] = [(host, TreeType.HOST)]
            else: # add leaves specifically for RL events on branches that lead to leaves
                if (parasite, TreeType.PARASITE) in temporal_graph:
                    print("3")
                    temporal_graph[(parasite, TreeType.PARASITE)].append((host, TreeType.HOST))
                else: # add leaves specifically for RL events on branches that lead to leaves
                    print("4")
                    temporal_graph[(parasite, TreeType.PARASITE)] = [(host, TreeType.HOST)]'''

            '''if (parasite, TreeType.PARASITE) in temporal_graph:
                if (host, TreeType.HOST) in temporal_graph:
                    if (parasite, TreeType.PARASITE) not in temporal_graph[(host, TreeType.HOST)]:
                        if parasite == 'g2':
                            print("g2")
                        temporal_graph[(parasite, TreeType.PARASITE)].append((host, TreeType.HOST))
                else:
                    temporal_graph[(parasite, TreeType.PARASITE)].append((host, TreeType.HOST))
            else: # add leaves specifically for RL events on branches that lead to leaves
                temporal_graph[(parasite, TreeType.PARASITE)] = [(host, TreeType.HOST)]'''
        # if the node_mapping is not a mapping onto the root of host tree, we add the second relation
        """if host_parent != 'h_root' and host_parent != 'p_root' and event_type != 'C' and event_type != 'L': # didn't exclude 'C' previously
            # Why are host_parents added ? -beta
            print("NOT many things: " + str(host_parent) + ", " + str(parasite))
            print(event_tuple)
            '''if (host, TreeType.HOST) not in temporal_graph[(parasite, TreeType.PARASITE)]: # BETA'''
            temporal_graph[(host_parent, TreeType.HOST)].append((parasite, TreeType.PARASITE))"""
        
        # if event is a transfer, then we add two more temporal relations
        if event_type == 'T' and add_strong_constraints:
            # get the mapping for the right child which is the transferred child
            right_child_mapping = event_tuple[2]
            right_child_parasite, right_child_host = right_child_mapping
            # since a transfer event is horizontal, we have these two implied relations
            temporal_graph[(parent[right_child_host], TreeType.HOST)].append((parasite, TreeType.PARASITE))

    for node_tuple in temporal_graph:
        # we need to make sure the associated value in the dictionary does not contain repeated node tuples
        temporal_graph[node_tuple] = uniquify(temporal_graph[node_tuple])
    return temporal_graph

# This is a topological sort based on depth-first-search
# https://en.wikipedia.org/wiki/Topological_sorting#Depth-first_search


def topological_order(temporal_graph):
    """
    :param temporal_graph: as described in the return type of build_temporal_graph
    :return: A dictionary in which a key is a node tuple (name, type) as described
        in build_temporal_graph and the value is a positive integer representing its topological ordering.
        The ordering numbers are consecutive values beginning at 1.
        If the graph has a cycle and the topological ordering therefore fails, this
        function returns None.
    """
    # the ordering of nodes starts at 1
    next_order = 1
    unvisited_nodes = OrderedDict.fromkeys(sorted(temporal_graph.keys()))
    # the visitng_nodes is used to detect cycles. If the visiting_nodes add an element that is already
    # in the list, then we have found a cycle
    visiting_nodes = set()
    ordering_dict = {}
    while unvisited_nodes:
        # removes the first node from unvisited_nodes
        start_node = unvisited_nodes.popitem(last=False)[0]
        has_cycle, next_order = topological_order_helper(start_node, next_order, visiting_nodes,
                               unvisited_nodes, temporal_graph, ordering_dict)
        if has_cycle: return None
    # reverse the ordering of the nodes
    for node_tuple in ordering_dict:
        ordering_dict[node_tuple] = next_order - ordering_dict[node_tuple]
    return ordering_dict


def topological_order_helper(start_node, start_order, visiting_nodes, unvisited_nodes, temporal_graph, ordering_dict):
    """
    :param start_node: is the starting node to explore the temporal_graph
    :param start_order: is the order we start to label the nodes with
    :param visiting_nodes: are nodes that are on the same path and are currently being explored
    :param unvisited_nodes: are nodes in temporal graph that have not been visited
    :param temporal graph: as described in the return type of build_temporal_graph
    :param ordering_dict: is the dictionary that contains labeled node tuples and their ordering as described
            in topological_order
    :return: a Boolean value that denotes whether the part of temporal graph reachable from start_node
             contains a cycle
    :return: the start order to be used by the remaing nodes of temporal graph that have not been labeled
    """
    next_order = start_order
    is_leaf = start_node not in temporal_graph or start_node[0].isdigit() # leaves should always be numbers for gene trees
    if is_leaf:
        return False, next_order
    else:
        has_cycle = start_node in visiting_nodes
        if has_cycle:
            return True, next_order
        visiting_nodes.add(start_node)
        child_nodes = sorted(temporal_graph[start_node])
        for child_node in child_nodes:
            # if the child_node is already labeled, we skip it
            if child_node in ordering_dict:
                continue
            if child_node in unvisited_nodes:
                unvisited_nodes.pop(child_node)
            has_cycle_child, next_order = topological_order_helper(child_node, next_order,  visiting_nodes,
                                            unvisited_nodes, temporal_graph, ordering_dict)
            # if we find a cycle, we stop the process
            if has_cycle_child: return True, next_order
        # if children are all labeled, we can label the start_node
        visiting_nodes.remove(start_node)
        ordering_dict[start_node] = next_order
        return False, next_order + 1


def populate_nodes_with_order(tree_node, tree_type, ordering_dict, leaf_order):
    """
    :param tree_node: the root node of the subtree we want to populate the temporal order information
    :param tree_type: the type of the tree
    :param ordering_dict: a dictionary that maps node tuples to their temporal order as described in topological_order
    :param leaf_order: the temporal order we should assign to the leaves of the tree
    """
    layout = NodeLayout()
    if tree_node.is_leaf():
        layout.col = leaf_order
        tree_node.layout = layout
    else:
        node_tuple = (tree_node.name, tree_type)
        layout.col = ordering_dict[node_tuple]
        tree_node.layout = layout
        populate_nodes_with_order(tree_node.left_node, tree_type, ordering_dict, leaf_order)
        populate_nodes_with_order(tree_node.right_node, tree_type, ordering_dict, leaf_order)


def _get_names_of_internal_nodes(tree):
    """
    :param: A host or parasite tree
    :return: A list of the names (strings) of the internal nodes in that tree
    """
    node_names = list()
    for edge_name in tree:
        edge_four_tuple = tree[edge_name]
        if not _is_leaf_edge(edge_four_tuple):
            node_names.append(_bottom_node(edge_four_tuple))
    return node_names


def _top_node(edge_four_tuple):
    """
    :param: 4-tuple of the form (top_vertex_name, bottom_vertex_name, child_edge1, child_edge2)
    :return: top_vertex_name
    """
    return edge_four_tuple[0]


def _bottom_node(edge_four_tuple):
    """
    :param: 4-tuple of the form (top_vertex_name, bottom_vertex_name, child_edge1, child_edge2)
    :return: bottom_vertex_name
    """
    return edge_four_tuple[1]


def _is_leaf_edge(edge_four_tuple):
    """
    :param: 4-tuple of the form (top_vertex_name, bottom_vertex_name, child_edge1, child_edge2)
    :return: True if child_edge1 = child_edge2 = None.
        This signifies that this edge terminates at a leaf.
    """
    return edge_four_tuple[3] is None
    
##################################################################
##################################################################
#################$ NECESSARY RECON FUNCTIONS #####################
##################################################################
##################################################################

"""
A Proposal for reconciliation graph and reconciliation interface.
"""

#__all__ = ['MappingNode', 'EventType', 'Event', 'Cospeciation', 'Duplication', 'Transfer', 'Loss', 'TipTip', 'Reconciliation', 'ReconGraph']

class MappingNode:
    """
    MappingNode is a node in the reconciliation graph that maps
    a parasite to a host.
    """
    def __init__(self, parasite_vertex: str, host_vertex: str):
        self._parasite = parasite_vertex
        self._host = host_vertex

    @property
    def parasite(self) -> str:
        return self._parasite

    @property
    def host(self) -> str:
        return self._host
    
    def __eq__(self, other):
        if type(self) == type(other):
            return self._parasite == other.parasite and self._host == other.host
        return False
    
    def __hash__(self):
        return hash((self._parasite, self._host))
    
    def __repr__(self):
        return "%s(%s, %s)" % (type(self).__name__, self._parasite, self._host)

class EventType(Enum):
    """
    EventType is the type of event node.
    """
    COSPECIATION = 1
    DUPLICATION = 2
    TRANSFER = 3
    LOSS = 4
    TIPTIP = 5
    ORIGIN = 6
    REARRANGEMENT = 7

class Event(ABC):
    """
    Base Class for events.
    """
    def __init__(self, freq: float = None):
        self._freq = freq

    @property
    def freq(self) -> float:
        return self._freq

    @property
    def event_type(self) -> EventType:
        """
        Returns the event type of this event.
        """
        raise NotImplementedError("Event is an abstract base class")

class TwoChildrenEvent(Event, ABC):
    """
    Base Class for events with two children.
    """
    def __init__(self, left: MappingNode, right: MappingNode, freq=None):
        """
        Creates a Cospeciation, Duplication, or Transfer event. The event node
        has two children MappingNodes called left and right.
        """
        super().__init__(freq)
        self._left: MappingNode = left
        self._right: MappingNode = right

    @property
    def left(self) -> MappingNode:
        return self._left

    @property
    def right(self) -> MappingNode:
        return self._right
    
    def __eq__(self, other):
        if type(self) == type(other):
            return self._left == other.left and self._right == other.right
        return False
    
    def __hash__(self):
        return hash((self._left, self._right))
    
    def __repr__(self):
        return "%s(%s, %s)" % (type(self).__name__, self._left, self._right)

class Cospeciation(TwoChildrenEvent):
    """
    Cospeciation is an event where two children of the parasite vertex
    maps to the two children of the host vertex. Cospeciation event has
    two child MappingNode named left and right.
    """

    @property
    def event_type(self) -> EventType:
        return EventType.COSPECIATION

class Duplication(TwoChildrenEvent):
    """
    Duplication is an event where two children of the parasite vertex
    maps to the same host vertex as their parent. Cospeciation event has
    two child MappingNode named left and right.
    """

    @property
    def event_type(self) -> EventType:
        return EventType.DUPLICATION

class Transfer(TwoChildrenEvent):
    """
    Transfer is an event where one child of the parasite vertex
    maps to the same host vertex as its parent, while the other child maps 
    to a host vertex not ancestrally related to the first one. Transfer event 
    has two child MappingNode named left and right. The right child is always
    the child that gets transferred.
    """

    @property
    def event_type(self) -> EventType:
        return EventType.TRANSFER

class Origin(TwoChildrenEvent):
    """
    Origin is an event where two children of the parasite vertex map to a single
    host vertex. There is no parent. Origin event has two child MappingNode 
    named left and right."""

    @property
    def event_type(self) -> EventType:
        return EventType.ORIGIN

class Loss(Event):
    """
    Loss is an event where the parasite vertex goes down to only one of 
    the lineage of the host vertex when the host vertex speciates. Loss event
    has one child MappingNode called child. The parasite of the child MappingNode
    (self.child.parasite) is the same as the parasite of the MappingNode that
    incurs the lost.
    """

    def __init__(self, child: MappingNode, freq=None):
        """
        Creates a Loss event. The lost event has one child MappingNode
        called child.
        """
        super().__init__(freq)
        self._child = child

    @property
    def child(self) -> MappingNode:
        return self._child
    
    def __eq__(self, other):
        return type(self) == type(other) and self._child == other.child
    
    def __hash__(self):
        return hash(self._child)
    
    def __repr__(self):
        return "%s(%s)" % (type(self).__name__, self._child)

    @property
    def event_type(self) -> EventType:
        return EventType.LOSS

# added by Michelle
class Rearrangement(Event):
    """
    Rearrangement is an event where the parasite vertex mapps to the same host
    vertex as the parent. 
    Rearrangement event has one child MappingNode called child. The parasite of 
    the child MappingNode (self.child.parasite) is the same as the parasite of 
    the MappingNode that incurs the rearrangement.
    """

    def __init__(self, child: MappingNode, freq=None):
        """
        Creates a Rearrangement event. The lost event has one child MappingNode
        called child.
        """
        super().__init__(freq)
        self._child = child

    @property
    def child(self) -> MappingNode:
        return self._child
    
    def __eq__(self, other):
        return type(self) == type(other) and self._child == other.child
    
    def __hash__(self):
        return hash(self._child)
    
    def __repr__(self):
        return "%s(%s)" % (type(self).__name__, self._child)

    @property
    def event_type(self) -> EventType:
        return EventType.REARRANGEMENT

class TipTip(Event):
    """
    TipTip (tip mapping) is an event where both the parasite and the host are 
    leaves. Biologically, this means the parasite and the host lives together
    at present time. The TipTip event is the sink of the reconciliation graph and
    has no children.
    """

    def __init__(self, freq: float = None):
        """
        Creates a TipTip (tip to tip mapping) event. The TipTip event has 
        no children mapping nodes.
        """
        super().__init__(freq)

    def __repr__(self):
        return "%s()" % type(self).__name__

    @property
    def event_type(self) -> EventType:
        return EventType.TIPTIP

class BranchReconciliation:
    """
    Reconciliation is a tree that represents an edge-to-edge mapping from 
    a parasite tree to a host tree. The Reconciliation starts with source
    MappingNode where the parasite first enters the host. Each MappingNode
    has one corresponding Event. Each event has up to two children MappingNode
    depending on its type. The leaves of the tree are the TipTip events which
    have no children MappingNode.
    """
    def __init__(self, source: MappingNode, initial_map: Dict[MappingNode, Event] = {}):
        """
        Creates a Reconciliation tree. The Reconciliation starts with source
        MappingNode where the parasite first enters the host.
        """
        self.source = source
        self._map = initial_map
        self._parasite_map: Dict[str, list] = {}
    
    def __eq__(self, other):
        return type(self) == type(other) and self._map == other._map
    
    def __repr__(self):
        return "%s(source=%s, %s)" % (type(self).__name__, self.source, self._map)
    
    def set_event(self, mapping: MappingNode, event: Event):
        """
        Set the event corresponding to mapping.
        """
        # if event.event_type is not EventType.LOSS:
        if mapping.parasite not in self._parasite_map:
            self._parasite_map[mapping.parasite] = [mapping]
        else:
            self._parasite_map[mapping.parasite].append(mapping)
        # add to map
        self._map[mapping] = event
    
    def event_of(self, mapping: MappingNode) -> Event:
        """
        Get the event corresponding to mapping.
        """
        return self._map[mapping]

    def mapping_of(self, parasite: str) -> MappingNode:
        """
        Returns the MappingNode of the parasite vertex. This is where that
        parasite vertex maps to. Never returns a mapping node that corresponds
        to a loss event.
        """
        return self._parasite_map[parasite]
    
    def is_sink(self, mapping: MappingNode) -> bool:
        """
        Whether the host and the parasite of the mapping are leaves of its
        phylogenic tree with a tip mapping.
        """
        return self._map[mapping].event_type is EventType.TIPTIP

    def __contains__(self, mapping: MappingNode):
        return True if mapping in self._map else False
    
    def save(self, path: str, metadata: Dict[str, str] = {}):
        """
        Save the Reconciliation to path, e.g. ``save('./reconname')``
        This will save the Reconciliation along with the metadata if metadata
        is specified, e.g. ``save('./reconname', {'created_at': 'Noon'})``.
        The metadata is only for book-keeping purposes and is not read when load.
        """
        raise NotImplementedError
    
    @staticmethod
    def load(path) -> 'Reconciliation':
        """
        Load a Reconciliation from path, e.g. ``Reconciliation.load('./reconname')``.
        """
        raise NotImplementedError

class Reconciliation:
    """
    Reconciliation is a tree that represents an edge-to-edge mapping from 
    a parasite tree to a host tree. The Reconciliation starts with source
    MappingNode where the parasite first enters the host. Each MappingNode
    has one corresponding Event. Each event has up to two children MappingNode
    depending on its type. The leaves of the tree are the TipTip events which
    have no children MappingNode.
    """
    def __init__(self, source: MappingNode, initial_map: Dict[MappingNode, Event] = {}):
        """
        Creates a Reconciliation tree. The Reconciliation starts with source
        MappingNode where the parasite first enters the host.
        """
        self.source = source
        self._map = initial_map
        self._parasite_map: Dict[str, MappingNode] = {}
    
    def __eq__(self, other):
        return type(self) == type(other) and self._map == other._map
    
    def __repr__(self):
        return "%s(source=%s, %s)" % (type(self).__name__, self.source, self._map)
    
    def set_event(self, mapping: MappingNode, event: Event):
        """
        Set the event corresponding to mapping.
        """
        if event.event_type is not EventType.LOSS:
            self._parasite_map[mapping.parasite] = mapping
        self._map[mapping] = event
    
    def event_of(self, mapping: MappingNode) -> Event:
        """
        Get the event corresponding to mapping.
        """
        return self._map[mapping]

    def mapping_of(self, parasite: str) -> MappingNode:
        """
        Returns the MappingNode of the parasite vertex. This is where that
        parasite vertex maps to. Never returns a mapping node that corresponds
        to a loss event.
        """
        return self._parasite_map[parasite]
    
    def is_sink(self, mapping: MappingNode) -> bool:
        """
        Whether the host and the parasite of the mapping are leaves of its
        phylogenic tree with a tip mapping.
        """
        return self._map[mapping].event_type is EventType.TIPTIP

    def __contains__(self, mapping: MappingNode):
        return True if mapping in self._map else False
    
    def save(self, path: str, metadata: Dict[str, str] = {}):
        """
        Save the Reconciliation to path, e.g. ``save('./reconname')``
        This will save the Reconciliation along with the metadata if metadata
        is specified, e.g. ``save('./reconname', {'created_at': 'Noon'})``.
        The metadata is only for book-keeping purposes and is not read when load.
        """
        raise NotImplementedError
    
    @staticmethod
    def load(path) -> 'Reconciliation':
        """
        Load a Reconciliation from path, e.g. ``Reconciliation.load('./reconname')``.
        """
        raise NotImplementedError

class ReconGraph:
    """
    ReconGraph is a directed acyclic graph that compactlyrepresent multiple 
    Reconciliations. The Reconciliation starts with a list of source MappingNode
    where the parasite first enters the host. Each MappingNode has a list of 
    corresponding Event. Each event has up to two children MappingNode
    depending on its type. The sinks of this graph are the TipTip events which
    have no children MappingNode.
    """
    def __init__(self, sources: List[MappingNode], initial_map: Dict[MappingNode, List[Event]] = {}):
        """
        Creates a ReconGraph. The ReconGraph starts with source
        MappingNode where the parasite first enters the host.
        """
        self.sources = sources
        self._map = initial_map
    
    def __eq__(self, other):
        return type(self) == type(other) and self._map == other._map
    
    def __repr__(self):
        return "%s(sources=%s, %s)" % (type(self).__name__, self.sources, self._map)
    
    def add_event(self, mapping: MappingNode, event: Event):
        """
        Add an event corresponding to mapping.
        """
        if mapping not in self._map:
            self._map[mapping] = []
        self._map[mapping].append(event)

    def set_events(self, mapping: MappingNode, events: List[Event]):
        """
        Set a mapping to correspond to a list of events.
        """
        self._map[mapping] = events
    
    def events_of(self, mapping: MappingNode) -> List[Event]:
        """
        Get a list of events corresponding to mapping.
        """
        return self._map[mapping].copy()
    
    def is_sink(self, mapping: MappingNode) -> bool:
        """
        Whether the host and the parasite of the mapping are leaves of its
        phylogenic tree with a tip mapping.
        """
        return len(self._map[mapping]) == 1 and self._map[mapping][0].event_type is EventType.TIPTIP
    
    def enumerate(self) -> List[Reconciliation]:
        """
        Returns all Reconciliation represented in this ReconGraph.
        """
        raise NotImplementedError

    def save(self, path: str, metadata: Dict[str, str] = {}):
        """
        Save the ReconGraph to path, e.g. ``save('./reconname')``
        This will save the reconciliation graph along with the metadata if metadata
        is specified, e.g. ``save('./reconname', {'created_at': 'Noon'})``.
        The metadata is only for book-keeping purposes and is not read when load.
        """
        raise NotImplementedError
    
    @staticmethod
    def load(path) -> 'ReconGraph':
        """
        Load a Reconciliation from path, e.g. ``Reconciliation.load('./reconname')``.
        """
        raise NotImplementedError

##################################################################
##################################################################
################### NECESSARY TREE FUNCTIONS #####################
##################################################################
##################################################################

"""
tree.py
Defines classes related to host and parasite nodes and trees
"""

class TreeType(Enum):
    """ Defines type of the tree """
    HOST = 1
    PARASITE = 2

class Track(Enum):
    HORIZONTAL = 1
    LOWER_VERTICAL = 2
    UPPER_VERTICAL = 3

class Node:
    """ Defines a node of a tree """
    def __init__(self, name):
        self.name = name            # String; name of this host or parasite
        self.left_node = None       # Node:  left child Node or None
        self.right_node = None      # Node:  right child Node or None
        self.parent_node = None     # Node:  parent Node or None
        self.layout = None          # NodeLayout object: layout of this node

    def is_leaf(self):
        """ returns True iff this node is a leaf/tip of the tree """
        return self.left_node is None and self.right_node is None

    def is_root(self):
        """ returns True iff this node is the root of the tree """
        return self.parent_node is None

    def __repr__(self):
        return str(self.name)

    def get_layout(self):
        """ returns the four values listed in NodeLayout"""
        layout = self.layout
        return layout.row, layout.col, layout.x, layout.y

    def set_layout(self, row=None, col=None, x=None, y=None):
        """Sets the layout"""
        layout = self.layout
        layout.row = row if row != None else layout.row
        layout.col = col if col != None else layout.col
        layout.x = x if x != None else layout.x
        layout.y = y if y != None else layout.y

    def get_and_update_track(self, track):
        """updates track number and returns previous track of host node"""
        if track == Track.HORIZONTAL:
            self.layout.h_track = self.layout.h_track + 1
            return self.layout.h_track - 1
        
        if track == Track.UPPER_VERTICAL:
            self.layout.upper_v_track += 1
            return self.layout.upper_v_track - 1
        
        if track == Track.LOWER_VERTICAL:
            self.layout.lower_v_track += 1
            return self.layout.lower_v_track - 1
    
    def update_count(self):
        self.layout.node_count += 1 


class NodeLayout:
    """ Defines node layout attributes for rendering """
    def __init__(self):
        self.row = None         # float: logical position of this Node in rendering

        # The self.col can be generated from a topological ordering of the temporal constraint graph
        self.col = None         # float: logical position of this Node in rendering
        self.x = None           # float: x-coordinate for rendering
        self.y = None           # float: y-coordinate for rendering
        self.upper_v_track = 1  # float: track number for upper vertical host edges
        self.lower_v_track = 1  # float: track number for lower vertical host edges
        self.h_track = 1        # int: track number for horizontal host edges
        self.node_count = 0     # int: Number of nodes mapped to node
        self.offset = 0         # int: Offset between tracks of a node


class Tree:
    """
    The Tree type defines a tree for use in rendering and other functions
    """
    def __init__(self):
        self.root_node = None       # Node:  Root Node of the Tree
        self.tree_type = None       # TreeType: HOST or PARASITE
        self.pos_dict = {}

    def leaf_list(self):
        """ Returns list of leaf Nodes from left to right. """
        return self._leaf_list_helper(self.root_node)

    def _leaf_list_helper(self, node):
        if node.is_leaf():
            return [node]
        list1 = self._leaf_list_helper(node.left_node)
        list2 = self._leaf_list_helper(node.right_node)
        list1.extend(list2)
        return list1

    def postorder_list(self):
        """ returns list of all Nodes in postorder """
        return self._postorder_list_helper(self.root_node)

    def name_to_node_dict(self):
        """
        Returns a dictionary whose keys are names (strings) and values are the
        nodes whose .name is that string.
        Use case:  A parasite p finds its corresponding host h and then uses this
        dictionary to get the h's node which contains, among other things, its layout.
        This allows the parasite p to set its layout based on that of the host.
        """
        ntn_dict = {}
        self._name_to_node_dict_helper(self.root_node, ntn_dict)
        return ntn_dict

    def _name_to_node_dict_helper(self, node, ntn_dict):
        ntn_dict[node.name] = node
        if node.is_leaf():
            return
        self._name_to_node_dict_helper(node.left_node, ntn_dict)
        self._name_to_node_dict_helper(node.right_node, ntn_dict)

    def _postorder_list_helper(self, node):
        if node.is_leaf():
            return [node]
        list1 = self._postorder_list_helper(node.left_node)
        list2 = self._postorder_list_helper(node.right_node)
        list1.extend(list2)
        list1.append(node)
        return list1

##################################################################
##################################################################
############ NECESSARY RENDER_SETTINGS FUNCTIONS #################
##################################################################
##################################################################

# py

VERTICAL_OFFSET = 0.3       # Offset for drawing parasite nodes above host nodes
COSPECIATION_OFFSET = .3    # Offest for drawing parasite nodes closer to host 
                            # nodes for speciation events
NODE_OFFSET = 0.3
TRACK_OFFSET = 0.3
TIP_TEXT_OFFSET_X = .3


# Colors
# Define new colors as 4-tuples of the form (r, g, b, 1) where
# r, g, b are values between 0 and 1 indicating the amount of red, green, and blue.
RED = (1, 0, 0, 1)


MAROON = (0.5, 0, 0, 1)
GREEN = (0, 0.5, 0, 1)

PURPLE = (0.5, 0, 0.5, 1)
BLACK = (0, 0, 0, 1)
GRAY = (0.5, 0.5, 0.5, 1)
WHITE = (1, 1, 1, 1)
PURPLE = (.843, .00, 1.0, 1)

BLUE = (.09, .216, .584, 1)
ROYAL_BLUE = (.3, .4, .9, 1)
CYAN = (.3, .9, .75, 1)
RED_BLUSH = (.882, .255, .412, 1)
PRETTY_YELLOW = (.882, .725, .255, 1)
ORANGE_ORANGE = (1.00, .502, 0, 1)

LEAF_NODE_COLOR = BLUE
COSPECIATION_NODE_COLOR = ORANGE_ORANGE
DUPLICATION_NODE_COLOR = PURPLE
TRANSFER_NODE_COLOR = RED_BLUSH
ORIGIN_NODE_COLOR = GREEN
REARRANGEMENT_NODE_COLOR = PRETTY_YELLOW
HOST_NODE_COLOR = BLACK
HOST_EDGE_COLOR = BLACK
PARASITE_EDGE_COLOR = ROYAL_BLUE
LOSS_EDGE_COLOR = GRAY

NORMAL_EDGE_COLOR = RED
LOSS_COLOR_SHAPE = RED, "x"

TRANSFER_TRANSPARENCY = 0.25

LEAF_NODE_SHAPE = "o"
COSPECIATION_NODE_SHAPE = "o"
DUPLICATION_NODE_SHAPE = "D"
TRANSFER_NODE_SHAPE = "s"
ORIGIN_NODE_SHAPE = "o"
REARRANGEMENT_NODE_SHAPE = "D"

TIP_ALIGNMENT = 'center'

CENTER_CONSTANT = 3 / 8

NODESIZE = 8
SPECIES_NODESIZE = 4
START_SIZE = -60
STEP_SIZE = 50
MIN_FONT_SIZE = 0
MAX_FONT_SIZE = .3
COUNT_OFFSET = 3
PUSHED_NODE_OFFSET = 0.5

INTERNAL_NODE_ALPHA = 0.7

HOST_NODE_BORDER_COLOR = WHITE
PARASITE_NODE_BORDER_COLOR = BLACK

TREE_TITLE = ""

LEGEND_ELEMENTS = [
       Line2D([0], [0], marker= ORIGIN_NODE_SHAPE, color='w', label = 'Origin',
              markerfacecolor=ORIGIN_NODE_COLOR, markersize=NODESIZE),
       Line2D([0], [0], marker= COSPECIATION_NODE_SHAPE, color='w', label='Cospeciation',
              markerfacecolor=COSPECIATION_NODE_COLOR, markersize=NODESIZE),
       Line2D([0], [0], marker=DUPLICATION_NODE_SHAPE, color='w', label='Duplication',
              markerfacecolor=DUPLICATION_NODE_COLOR, markersize=NODESIZE),
       Line2D([0], [0], marker=TRANSFER_NODE_SHAPE, color='w', label='Transfer',
              markerfacecolor=TRANSFER_NODE_COLOR, markersize=NODESIZE),
       Line2D([0], [0], marker=REARRANGEMENT_NODE_SHAPE, color='w', label='Rearrangement',
              markerfacecolor=REARRANGEMENT_NODE_COLOR, markersize=NODESIZE),
       LineCollection([[(0, 0)]], linestyles=['dashed'],
              colors=[LOSS_EDGE_COLOR], label='Loss')
]

UP_ARROW_ROTATION = 0
DOWN_ARROW_ROTATION = 180

##################################################################
##################################################################
############# NECESSARY PLOT_TOOLS FUNCTIONS #####################
##################################################################
##################################################################

"""
py
Plotting tools using matplotlib
"""

# If matplotlib doesn't pop up a window, force it to use tkinter backend
# matplotlib.use("tkagg")
LINEWIDTH = 1
TEXTWIDTH = .3
BORDER_WIDTH = 1.2

LINE_Z_ORDER = 0
DOT_Z_ORDER = 1
TEXT_Z_ORDER = 2

SIZE = 6
TRANSFERSIZE = 10

FONTSIZE = 12

DEFAULT_VERTICAL_ALIGNMENT = 'bottom'
DEFAULT_VERTICAL_ALIGNMENT_2 = 'top'
DEFAULT_HORIZONTAL_ALIGNMENT = 'right'
DEFAULT_LOCATION = 'lower left' #'upper left'
DEFAULT_LINESTYLE = '-'
DEFAULT_TRIANGLE_LINESTYLE = 'None'
DEFAULT_DOT_MARKER = 'o'
CENTER = 'center'

def transparent_color(col: Tuple[int, int, int, float], alpha: float):
    return col[0:3] + (alpha,)

class Position(NamedTuple):
    x: int
    y: int

class FigureWrapper:
    """ Class definining plotting methods """
    def __init__(self, title: str, axes: Union[plt.Axes, None] = None):
        """
        If axes is specified, draw on axes instead.
        """
        if axes is None:
            self.fig = plt.figure()
            self.axis = self.fig.subplots(1, 1) # creates a figure with one Axes (plot)
        else:
            self.fig = axes.get_figure()
            self.axis = axes
        self.axis.autoscale()
        self.axis.margins(0.1)
        self.axis.axis("off")
        self.axis.set_title(title)

    def set_legend(self, legend_elements: list, loc: str = DEFAULT_LOCATION, fontsize: int = FONTSIZE, title: str = None):
        """
        create legend
        """
        self.axis.legend(handles=legend_elements, loc=loc, title=title, fontsize=fontsize, title_fontsize=fontsize)

    def line(self, point_1: Position, point_2: Position, col: tuple = BLACK, linestyle: str = DEFAULT_LINESTYLE):
        """
        Draw line from point p1 to p2
        """
        x_1, y_1 = point_1
        x_2, y_2 = point_2
        self.axis.plot([x_1, x_2], [y_1, y_2], color=col, linewidth=LINEWIDTH, linestyle=linestyle, zorder=LINE_Z_ORDER)

    def dot(self, point: Position, marker: str = DEFAULT_DOT_MARKER, col: tuple = BLACK):
        """
        Plot dot at point p
        """
        self.axis.plot(point.x, point.y, marker, color=col, zorder=DOT_Z_ORDER)

    def halfDot(self, point: Position, marker: str = DEFAULT_DOT_MARKER, col: tuple = BLACK):
        """
        Plot slightly smaller dot at point p
        """
        self.axis.plot(point.x, point.y, marker, color=col, zorder = DOT_Z_ORDER, markersize=SPECIES_NODESIZE)

    def text(self, point: tuple, string: str, col: tuple = RED, size=FONTSIZE, h_a: str = DEFAULT_HORIZONTAL_ALIGNMENT):
        x, y = point
        self.axis.text(x, y, string, color=col, fontsize = size, horizontalalignment=h_a, verticalalignment=DEFAULT_VERTICAL_ALIGNMENT_2)

    def text_v2(self, point: tuple, text: str, col: tuple = BLACK, size: float = SIZE, vertical_alignment: str = DEFAULT_VERTICAL_ALIGNMENT, border_col: tuple = None):
        """
        Plot text string s at point p in monospace font
        """
        if text is not None:
            if vertical_alignment == CENTER:
                point = (point[0], point[1] - size * CENTER_CONSTANT)

            mono_property = font_manager.FontProperties(family='monospace')
            tp = TextPath(point, text, size=size, prop=mono_property)
            path_patch = PathPatch(tp, color=col, linewidth = TEXTWIDTH, zorder=TEXT_Z_ORDER)
            if border_col:
                path_patch.set_path_effects([PathEffects.withStroke(linewidth=BORDER_WIDTH, foreground=border_col)])
            self.fig.gca().add_patch(path_patch)
    
    def triangle(self, point: Position, col: tuple = BLACK, markersize: int = TRANSFERSIZE, rotation: float = UP_ARROW_ROTATION):
        """
        Draws a triangle in the desired position
        """
        self.axis.plot(point.x, point.y, color=col, marker=(3, 0, rotation), markersize=TRANSFERSIZE, linestyle=DEFAULT_TRIANGLE_LINESTYLE)

    def arrow_segment(self, point_1: Position, point_2: Position, col: tuple = BLACK):
        """
        Draws a line from point 1 to point 2 with an arrow in the middle
        """
        plt.plot((point_1.x, point_2.x), (point_1.y, point_2.y), linewidth=2, color=col)
        arrow_length_x, arrow_length_y = (point_2.x - point_1.x) / 2, (point_2.y - point_1.y) / 2
        plt.arrow(point_1.x, point_1.y, arrow_length_x, arrow_length_y, linewidth=2,
                  head_width=0.3, head_length=0.5, facecolor=col, edgecolor=col,
                  length_includes_head=False)


    def show(self):
        """
        Display figure
        """
        plt.figure(self.fig.number)
        plt.show()

    def save(self, filename: str):
        """
        Save figure to file
        """
        self.fig.savefig(filename)


##################################################################
##################################################################
########### NECESSARY XENOGI TRANSLATION FUNCTIONS ###############
##################################################################
##################################################################

from xenoGI.xenoGI import genomes, parameters, loadTreeRelatedData, loadGenomeRelatedData, families
from xenoGI.Family import *

#-----------------------------------------------------------------
# OriginFamiy dtlorMprD -> Reconstruction Dictionary
#-----------------------------------------------------------------

def addAppend(key, value, dictionary):
    """ adds or appends a value to a dictionary, depending on if the key exists
    """
    # if key is in dictionary, append
    if key in dictionary:
        dictionary[key].append(value)
    # if key is not in dictionary, we add it.
    else:
        dictionary[key] = [value]

def addOneKey(genesO:genomes.genes, geneTree:Tree, MprD:dict, mprdKey:tuple, outDict:dict, geneToKeyD:dict):
    """ adds a single entry to the reconciliation dictionary (outDict)
    :param genesO:  genesO object, called 'genesO' in interactive analysis
    :param geneTree:  Tree object representing the gene tree
    :param MprD:  most probable reconciliation dictionary
    :param mprdKey:  tuple representing the key for the MprD that we should be working with
    :param outDict:  reconciliation dictionary being edited
    :param geneToKeyD:  dictionary mapping gene names to the reconciliation keys
    :return geneToKeyD (updated)
    """
    # only run if the key we're looking for is in the dictionary
    if mprdKey in MprD:
        # unpack mprd entry
        geneLocation,geneNB = mprdKey
        reconKey = ()
        eventsOrdered = MprD[mprdKey][::-1]
        # check all events
        for eventT in eventsOrdered:
            event,speciesLocation,stNB,locus = eventT
            # create the mprdKey
            reconKey = (geneLocation, speciesLocation)
            reconValue = ()
            # L and R events have only one child (the next event)
            if event == 'L' or event == 'R':
                child = geneToKeyD[geneLocation]
                reconValue = (event, child, None)
            # DTO events have two children
            else:
                # check for children
                rightChild, leftChild = geneTree.children(geneLocation)
                leftChildKey = geneToKeyD[leftChild]
                rightChildKey = geneToKeyD[rightChild]
                reconValue = (event, leftChildKey, rightChildKey)
            # add event, and key
            addAppend(reconKey, reconValue, outDict)
            geneToKeyD[geneLocation] = reconKey

def mprdToReconHelper(genesO:genomes.genes, geneTree:Tree, MprD:dict, node:str, outDict:dict):
    """ creates the reconcilliation dictionary by recursing through enries in MprD.
        returns geneToOutKey dictionary, because higher levels in the tree are 
        dependent on this dictionary: to figure out what proper keys are for
        higher-levels, we need the keys from the lower levels, (including leaves)
    :param genesO:  genesO object, called 'genesO' in interactive analysis
    :param geneTree:  Tree object representing the gene tree
    :param MprD:  most probable reconciliation dictionary
    :param node: str representing the node to start at and recurse down from
    :param outDict: reconciliation dictionary being edited
    :return geneToOutKey, a dictionary mapping genes to their keys. This is 
            required for consistent naming
    """
    # at a tip
    key = ()
    value = ()
    # given a gene (as in, a child from geneTreeO.children(node)), give the associated key in outKey
    geneToOutKey = {}
    if geneTree.isLeaf(node): # give none, only one option
        # create the key for the leaf
        species = genesO.numToStrainName(int(node))
        key = (node, species)
        # add the key to dictionary mapping genes to appropriate keys
        geneToOutKey[node] = key
        if (node, 'b') in MprD: # may be a loss or rearrangement
            addOneKey(genesO, geneTree, MprD, (node, 'b'), outDict, geneToOutKey)
        value = ('C', (None, None), (None, None))
        addAppend(key, value, outDict)
    else: # go through both children
        leftChild, rightChild = geneTree.children(node)
        # add the appropriate out keys to the dictionary, and recurse
        geneToOutKey = {**geneToOutKey, **mprdToReconHelper(genesO, geneTree, MprD, leftChild, outDict)}
        geneToOutKey = {**geneToOutKey, **mprdToReconHelper(genesO, geneTree, MprD, rightChild, outDict)}
        # check for n
        mprdKey1 = (node, 'n')
        addOneKey(genesO, geneTree, MprD, mprdKey1, outDict, geneToOutKey)
        # check for b
        mprdKey2 = (node, 'b')
        addOneKey(genesO, geneTree, MprD, mprdKey2, outDict, geneToOutKey)
    return geneToOutKey

def mprdToRecon(genesO:genomes.genes, geneTree:Tree, MprD:dict):
    """ Translate an MprD (most probable recon dictionary) to a reconciliation 
        dictionary format that render takes in. 
    :param genesO:  genesO object, called 'genesO' in interactive analysis
    :param geneTree:  Rtree object representing the gene tree
    :param MprD:  the most probable reconciliation dictionary (ofam.dtlorMprD)
    :return reconciliation dictionary
    """
    # get the root
    root = geneTree.rootNode
    # initialize empty dictionary
    dictionary = {}
    # edit the dictionary
    mprdToReconHelper(genesO, geneTree, MprD, root, dictionary)
    return dictionary

#-----------------------------------------------------------------
# Tree things -> Better tree things
#-----------------------------------------------------------------

def getSpecies(reconDict:dict):
    """ Pulls the second value from each key in a reconstruction dictionary, yielding
        all species involved in the reconstruction.
    :param reconDict:  reconstruction dictionary
    :return List of species mentioned in the reconstruction
    """
    # initialize list
    species = []
    for key in reconDict:
        geneLocation, speciesLocation = key
        species.append(speciesLocation)
    return species

def catchLeafDictionary(dictionary:dict, replacementValue:str):
    """ Host(species)/Parasite(gene) trees containing only a single entry have
        a root parent's name that is not h_root or p_root, which causes errors
        farther down the line. This function catches those dictionaries, and sets them
        to a non-erroring state. (changes root parent to replacementValue)
    :param dictionary:  dictionary created by calling "Rtree.createDtlorD" on a given Rtree
    :param replacementValue:  either 'h_root' for hosts or 'p_root' for parasites
    :return nothing
    """
    # if dictionary is a leaf
    if len(dictionary) == 1:
        # should only be one key
        onlyKey = list(dictionary.keys())[0]
        # set the 'parent' of the only entry to replacementValue
        parent, node, leftChild, rightChild = dictionary[onlyKey]
        dictionary[onlyKey] = (replacementValue, node, leftChild, rightChild)

#-----------------------------------------------------------------
# Origin family -> figure
#-----------------------------------------------------------------

def ofamRender(genesO:genomes.genes, originFamilyObject:Family.Families, speciesTreeObject:Tree.Rtree, number:int, showLegend: bool = False, showInternalLabels: bool = True):
    """ Renders a reconciliation using 'render', given an origin family number
    :param genesO:  genesO object, called 'genesO' in interactive analysis
    :param originFamilyObject:  Family object containing all origin families, called 'originFamiliesO' in interactive analysis
    :param speciesTreeObject:  Host / species tree. Called 'speciesRtreeO' in interactive analysis
    :param number:  integer number for the origin family that reconstruction should be generated on.
    :return Figure Object
    """
    # get the origin family, and get its reconstruction dictionary
    ofam = originFamilyObject.getFamily(number)
    mprD = ofam.dtlorMprD
    reconDict = mprdToRecon(genesO, ofam.geneTreeO, mprD)
    
    # get the smallest possible tree that the reconstruction can be mapped on
    #   saves space by only using the portion of the species graph required
    speciesList = getSpecies(reconDict)
    mrca = speciesTreeObject.findMrca(speciesList)
    subtree = speciesTreeObject.subtree(mrca)
    
    # change trees into dictionaries, fixing them if we are working with a single species
    hostDict = subtree.createDtlorD(True)
    catchLeafDictionary(hostDict, 'h_root')
    parasiteDict = ofam.geneTreeO.createDtlorD(False)
    catchLeafDictionary(parasiteDict, 'p_root')

    # render and return the figure
    outputFigure = render(hostDict, parasiteDict, reconDict, show_legend = showLegend, show_internal_labels= showInternalLabels)
    return outputFigure

def plotOriginFamily(genesO, number:int):
    """ Wrapper function for ofamRender, should work in interactiveAnalysis
    :param number:  integer number for the origin family that reconstruction should be generated on.
    :return Figure Object
    """
    # render and return
    figure = ofamRender(genesO, originFamiliesO, speciesRtreeO, number)
    return figure

def setup():
    # creates genesO, speciesRtreeO, and originFamiliesO for plotOriginFamily
    # STOLEN FROM INTERACTIVE ANALYSIS
    paramFN = 'params.py'
    global paramD
    paramD = parameters.createParametersD(parameters.baseParamStr,paramFN)

    global genesO
    global speciesRtreeO
    global originFamiliesO
    strainNamesT,genesO,geneOrderD = loadGenomeRelatedData(paramD)
    speciesRtreeO,subtreeD = loadTreeRelatedData(paramD['speciesTreeFN'])
    genesO.initializeGeneInfoD(paramD['geneInfoFN'],strainNamesT)
    genesO.initializeGeneNumToNameD(paramD['geneInfoFN'],strainNamesT)
    originFamiliesO = families.readFamilies(paramD['originFamilyFN'],speciesRtreeO,genesO,"origin")

