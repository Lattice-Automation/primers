# primers

```python
from primers import pcr, Primer

# create a FWD and REV primer
ps = pcr("ATGGATGGTAGAGATAGATGG",
    add_fwd="GGTAGGTAGAT", add_rev="GGTTTTAGGATAGAT", add_min=10, add_max=25,
    opt_tm=55.0, opt_len=30)

ps[0] # Primer("GGTAGGTAGATATGGATGGT", tm=62.3)
ps[1] # Primer("GGTTTTAGGATAGATCCTACTATCT", tm=64.0)
```
