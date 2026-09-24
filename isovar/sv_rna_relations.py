"""How a reconstructed RNA junction relates to a nominated SV adjacency."""

# Linkage of a path to the nominated adjacency, strongest first. The first
# three are event-linked; the last two may arise without the rearrangement.
LINKAGE_STATUSES = (
    "breakpoint_junction",  # The RNA join reproduces the DNA adjacency.
    "event_compatible_junction",  # Donor-side to acceptor-side join, e.g. spliced.
    "breakpoint_clip_partner_unplaced",  # Unaligned sequence at a breakpoint.
    # Crosses the event, but ordinary forward splicing (or read-through) of
    # the unrearranged reference could make it: any such non-breakpoint join,
    # and a breakpoint join between annotated splice sites.
    "splice_ambiguous_event_junction",
    "regional_novel_junction",  # Unannotated join not crossing the event.
)

# Relations which rest on the rearrangement (the "event-linked" statuses).
EVENT_LINKED_RELATIONS = LINKAGE_STATUSES[:3]

# Relations which cross the event, including splice-ambiguous joins: the
# junctions an exploratory ORF or event-linked comparison must cross.
EVENT_CROSSING_RELATIONS = LINKAGE_STATUSES[:4]
