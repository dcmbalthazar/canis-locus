def interpret_phenotype(geno):
    alerts = []

    # 1) Locus E — vermelho australiano (recessive red)
    if geno["E"].count("e") == 2:
        phenotype = "vermelho australiano"
        if geno["M"].count("M") == 1:
            alerts.append("Porta merle não visível")
    else:
        # 2) Locus A — sable / tricolor / sólido
        A = geno["A"]
        if "Ay" in A:
            pattern = "sable"
        elif "at" in A:
            pattern = "tricolor"
        else:
            pattern = None  # sólido

        # 3) Locus B — preto / marrom
        base = "marrom" if geno["B"].count("b") == 2 else "preto"

        # 4) Locus D — diluição
        if geno["D"].count("d") == 2:
            if base == "preto":
                base = "azul"
            elif base == "marrom":
                base = "lilac"

        # 5) Locus M — merle
        is_merle = geno["M"].count("M") == 1

        # 6) Montar nome CBKC-style (uso no Brasil)
        if pattern == "sable":
            phenotype = "sable merle" if is_merle else "sable"
        elif pattern == "tricolor":
            if is_merle:
                # aqui, "preto merle" vira "azul merle"
                if base == "preto":
                    phenotype = "azul merle tricolor"
                else:
                    phenotype = f"{base} merle tricolor"
            else:
                phenotype = f"{base} tricolor"
        else:
            if is_merle:
                # no uso comum, "preto merle" é chamado de "azul merle"
                if base == "preto":
                    phenotype = "azul merle"
                else:
                    phenotype = f"{base} merle"
            else:
                phenotype = base

    # 7) Locus S — branco
    sp = geno["S"].count("sp")
    if sp == 2:
        phenotype += " (branco excessivo)"
        alerts.append("Risco aumentado de surdez/alterações oculares")
    elif sp == 1:
        phenotype += " (branco moderado)"
    else:
        phenotype += " (branco mínimo)"

    return phenotype.strip(), alerts
