class OwnableObject:
    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        kwargs["user"] = self.request.user
        return kwargs

    def attach_owner(self, form):
        self.object = form.save(commit=False)
        if not self.object.owner:
            self.object.owner = self.request.user
        self.object.save()
