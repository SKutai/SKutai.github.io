import { CommonModule } from '@angular/common';
import { Component } from '@angular/core';

@Component({
  selector: 'app-contact',
  standalone: true,
  imports: [CommonModule],
  templateUrl: './contact.component.html',
  styleUrls: ['./contact.component.scss'],
})
export class ContactComponent {
  focusOnEmail: boolean = false;

  constructor() {}

  ngOnInit() {}

  onFocus() {
    this.focusOnEmail = true;
  }

  offFocus() {
    this.focusOnEmail = false;
  }

}
