from tabulate import tabulate

def tabulate_de_results(X, m1, m2):
    # Step 1: Get all inner keys across all X[i]

    #all_inner_keys = sorted({k for subdict in X.values() for k in subdict})

    def custom_key(k):
        if k.endswith("+"):
            return (1, int(k[:-1]))  # "+"-keys go after, sort by numeric part
        else:
            return (0, int(k))       # regular keys first, numeric sort

    all_inner_keys = sorted({k for subdict in X.values() for k in subdict}, key=custom_key)




    # Step 2: Build the table
    table = []
    for outer_key in X:
        row = [outer_key]
        for inner_key in all_inner_keys:
            val = X[outer_key].get(inner_key)
            if val is not None and len(val) > max(m1, m2):
                row.append((val[m1], val[m2]))
            else:
                row.append(None)
        table.append(row)

    # Step 3: Sort rows by number of filled entries (tuples ≠ None)
    table.sort(key=lambda row: sum(v is not None for v in row[1:]), reverse=True)

    # Step 4: Print using tabulate (stringify tuples)
    headers = ["Key"] + all_inner_keys
    formatted_table = [
        [row[0]] + [f"({v[0]:.3f}, {v[1]:.3f})" if isinstance(v, tuple) else "" for v in row[1:]]
        for row in table
    ]
    print(tabulate(formatted_table, headers=headers, tablefmt="grid"))
