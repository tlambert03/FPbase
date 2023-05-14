def long_blurb(self, withbright=False, withbleach=False):
    blurb = self.name
    switch = self.get_switch_type_display().lower()
    if switch == "basic":
        switch += " (constitutively fluorescent)"
    blurb += " is a{}{}{} fluorescent protein".format(
        "n" if switch.startswith(("a", "e", "i", "o", "u")) else "",
        " " + switch if switch else "",
        " " + self.color.lower() if self.color else "",
    )
    if self.primary_reference:
        blurb += f" published in {self.primary_reference.year}"
    if self.parent_organism:
        blurb += f", derived from {self.parent_organism}. "
    else:
        blurb += ".  "

    if self.default_state:
        bright = None
        if withbright:
            if self.default_state.local_brightness > 1.5:
                bright = "much brighter than"
            elif self.default_state.local_brightness > 1.1:
                bright = "slightly brighter than"
            elif self.default_state.local_brightness > 0.9:
                bright = "of average brightness compared to"
            elif self.default_state.local_brightness >= 0.5:
                bright = "slightly dimmer than"
            elif self.default_state.local_brightness < 0.5:
                bright = "much dimmer than"
        if bright:
            blurb += f"It is {bright} other proteins in the database with similar emission spectra"

        bleach = None
        if withbleach:
            BM = self.default_state.bleach_measurements.first()
            if BM:
                if BM.rate > 300:
                    bleach = "excellent"
                elif BM.rate > 200:
                    bleach = "very good"
                elif BM.rate > 100:
                    bleach = "decent"
                elif BM.rate > 50:
                    bleach = "relatively poor"
                elif BM.rate < 50:
                    bleach = "very poor"
        if bleach:
            if bright:
                if (BM.rate > 100 and self.default_state.local_brightness < 0.9) or (
                    BM.rate <= 100 and self.default_state.local_brightness > 0.9
                ):
                    join = " but"
                else:
                    join = " and"
            else:
                join = "It"
            blurb += f"{join} has {bleach} photostability.  "
        elif bright or bleach:
            blurb += ".  "

        mature = None
        M = self.default_state.maturation
        if M:
            if M <= 15:
                mature = "very rapidly-maturing"
            elif M < 50:
                mature = "rapidly-maturing"
            elif M < 90:
                mature = "somewhat slowly-maturing"
            elif M >= 90:
                mature = ""
        if mature:
            blurb += "It is reported to be a {} {}".format(mature, self.get_agg_display().lower() or "protein")

        acid = None
        A = self.default_state.pka
        if A:
            if A > 6:
                acid = "high"
            elif A > 5:
                acid = "moderate"
            elif A >= 4:
                acid = "low"
            elif A < 4:
                acid = "very low"
        if acid:
            if mature:
                blurb += " with "
            else:
                blurb += "It has "
            blurb += f"{acid} acid sensitivity."
        if mature and not (acid or A):
            blurb += "."
    else:
        if self.get_agg_display():
            blurb += f"It is reported to be a {self.get_agg_display().lower()}."
    if self.cofactor:
        blurb += f" It requires the cofactor {self.get_cofactor_display().lower()} for fluorescence."

    return blurb.strip()
